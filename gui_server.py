from flask import Flask, render_template_string, request, jsonify, send_from_directory
import subprocess
import json
import os
import sys
import threading
import time
from pathlib import Path

app = Flask(__name__)

# Global variables to track simulation progress
simulation_running = False
simulation_progress = 0
simulation_status = "Idle"
simulation_log = []

def run_simulation_with_progress(params):
    """Run the BEC simulation with progress tracking"""
    global simulation_running, simulation_progress, simulation_status, simulation_log
    
    simulation_running = True
    simulation_progress = 0
    simulation_status = "Initializing"
    simulation_log = ["Starting simulation..."]
    
    try:
        # Change to the src directory
        src_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src')
        os.chdir(src_dir)
        
        # Activate conda environment and run the simulation
        # First, try to run the simulation using the current Python environment
        # In a real scenario, you might need to activate a specific conda environment
        cmd = [sys.executable, 'main.py']
        
        # Add parameters to command
        if 'atomic_parameters' in params:
            atomic = params['atomic_parameters']
            cmd.extend([
                '--relative_atomic_mass', str(atomic.get('relative_atomic_mass', 87)),
                '--scattering_length', str(atomic.get('scattering_length', -10.0)),
                '--atom_number', str(atomic.get('atom_number', 2500))
            ])
        
        if 'time_structure' in params:
            time_struct = params['time_structure']
            cmd.extend([
                '--duration', str(time_struct.get('duration', 200.0)),
                '--dt', str(time_struct.get('dt', 0.5)),
                '--sampling_interval', str(time_struct.get('sampling_interval', 10))
            ])
        
        if 'grid_structure' in params:
            grid = params['grid_structure']
            cmd.extend([
                '--radius_xy', str(grid['radius_xy'][0]), str(grid['radius_xy'][1]),
                '--number_xy', str(grid['number_xy'][0]), str(grid['number_xy'][1])
            ])
        
        if 'potential_trap' in params:
            trap = params['potential_trap']
            cmd.extend([
                '--omega_trap', str(trap.get('omega_trap', 0.031415926535897934)),
                '--center_trap', str(trap['center_trap'][0]), str(trap['center_trap'][1]),
                '--beta', str(trap.get('beta', -0.1)),
                '--r_0', str(trap.get('r_0', 30.0)),
                '--omega_trap_z', str(trap.get('omega_trap_z', 0.015707963267948967))
            ])
        
        if 'initial_wavepacket' in params:
            wavepacket = params['initial_wavepacket']
            cmd.extend([
                '--center_bec', str(wavepacket['center_bec'][0]), str(wavepacket['center_bec'][1]),
                '--omega_bec', str(wavepacket.get('omega_bec', 0.031415926535897934)),
                '--velocity', str(wavepacket['velocity'][0]), str(wavepacket['velocity'][1]),
                '--angular_momentum_bec', str(wavepacket['angular_momentum_bec'][0]), str(wavepacket['angular_momentum_bec'][1])
            ])
        
        if 'simulation_settings' in params:
            settings = params['simulation_settings']
            if settings.get('imaginary_time', False):
                cmd.append('--imaginary_time')
            if settings.get('video', False):
                cmd.append('--video')
            if settings.get('mechanics', False):
                cmd.append('--mechanics')
            if settings.get('figure', False):
                cmd.append('--figure')
        
        # Log the command being executed
        simulation_log.append(f"Executing command: {' '.join(cmd)}")
        
        # Execute the simulation
        simulation_progress = 5
        simulation_status = "Activating environment and starting simulation..."
        
        # Use Popen to run the command and capture output
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            bufsize=1
        )
        
        # Read output in real-time
        for line in iter(process.stdout.readline, ''):
            simulation_log.append(line.strip())
            # Update progress based on output from the simulation
            line_lower = line.lower()
            if "parameters:" in line_lower:
                simulation_progress = 10
                simulation_status = "Processing parameters..."
            elif "time evolution loop" in line_lower or "time evolution" in line_lower:
                simulation_progress = 30
                simulation_status = "Running time evolution..."
            elif "sampling" in line_lower or "sample" in line_lower:
                simulation_progress = 50
                simulation_status = "Sampling wavefunction..."
            elif "saving figures" in line_lower or "saving figures" in line_lower:
                simulation_progress = 80
                simulation_status = "Saving results..."
            elif "done!" in line_lower:
                simulation_progress = 100
                simulation_status = "Completed"
        
        # Wait for process to complete
        process.wait()
        
        simulation_status = "Completed" if process.returncode == 0 else f"Failed with code {process.returncode}"
        
    except Exception as e:
        simulation_status = f"Error: {str(e)}"
        simulation_log.append(f"Error: {str(e)}")
    finally:
        simulation_running = False
        if simulation_status != "Completed":
            simulation_progress = 100


@app.route('/')
def index():
    """Serve the main GUI page"""
    gui_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gui.html')
    with open(gui_path, 'r', encoding='utf-8') as f:
        html_content = f.read()
    return render_template_string(html_content)


@app.route('/api/defaults', methods=['GET'])
def get_defaults():
    """Get default parameters"""
    defaults_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src', 'defaults.json')
    with open(defaults_path, 'r') as f:
        defaults = json.load(f)
    return jsonify(defaults)


@app.route('/api/defaults', methods=['POST'])
def save_defaults():
    """Save default parameters"""
    data = request.json
    defaults_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src', 'defaults.json')
    with open(defaults_path, 'w') as f:
        json.dump(data, f, indent=2)
    return jsonify({"status": "success"})


@app.route('/api/simulation/start', methods=['POST'])
def start_simulation():
    """Start the simulation"""
    global simulation_running
    
    if simulation_running:
        return jsonify({"status": "error", "message": "Simulation already running"})
    
    params = request.json
    
    # Start the simulation in a separate thread
    thread = threading.Thread(target=run_simulation_with_progress, args=(params,))
    thread.daemon = True
    thread.start()
    
    return jsonify({"status": "started"})


@app.route('/api/simulation/status', methods=['GET'])
def get_simulation_status():
    """Get current simulation status"""
    global simulation_running, simulation_progress, simulation_status, simulation_log
    return jsonify({
        "running": simulation_running,
        "progress": simulation_progress,
        "status": simulation_status,
        "log": simulation_log[-20:] if simulation_log else []  # Return last 20 log entries
    })


@app.route('/api/simulation/stop', methods=['POST'])
def stop_simulation():
    """Stop the simulation"""
    global simulation_running, simulation_status
    simulation_running = False
    simulation_status = "Stopped by user"
    return jsonify({"status": "stopped"})


@app.route('/api/results/images', methods=['GET'])
def get_result_images():
    """Get available result images"""
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src', 'output')
    if not os.path.exists(output_dir):
        return jsonify({"images": [], "categorized": {}})
    
    image_extensions = ['.png', '.jpg', '.jpeg', '.gif']
    images = []
    
    for filename in os.listdir(output_dir):
        if any(filename.lower().endswith(ext) for ext in image_extensions):
            images.append(filename)
    
    # Categorize images
    categorized = {
        'density': [img for img in images if 'density' in img.lower()],
        'phase': [img for img in images if 'phase' in img.lower()],
        'flow': [img for img in images if 'flow' in img.lower() or 'field' in img.lower()],
        'real': [img for img in images if 'real' in img.lower()],
        'quantity': [img for img in images if 'quantity' in img.lower() or 
                     any(q in img.lower() for q in ['Iz', 'Lz', 'omega', 'ang'])],
        'other': []
    }
    
    # Separate final state images for easier access
    categorized['density_final'] = [img for img in categorized['density'] if 'final' in img.lower()]
    categorized['phase_final'] = [img for img in categorized['phase'] if 'final' in img.lower()]
    categorized['flow_final'] = [img for img in categorized['flow'] if 'final' in img.lower()]
    categorized['real_final'] = [img for img in categorized['real'] if 'final' in img.lower()]
    
    # Add any uncategorized images to 'other'
    all_categorized = []
    for category in categorized.values():
        all_categorized.extend(category)
    
    categorized['other'] = [img for img in images if img not in all_categorized]
    
    return jsonify({"images": images, "categorized": categorized})


@app.route('/results/<filename>')
def send_result_image(filename):
    """Serve result images"""
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src', 'output')
    return send_from_directory(output_dir, filename)


if __name__ == '__main__':
    print("Starting BEC Dynamics 2D GUI server...")
    print("Open your browser and go to http://127.0.0.1:5000")
    app.run(debug=True, host='127.0.0.1', port=5000)