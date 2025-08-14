import tkinter as tk
from tkinter import messagebox

class GUI:
    def __init__(self, master, language):
        self.master = master
        self.language = language
        self.master.title(self.language['title'])

        self.label = tk.Label(master, text=self.language['welcome_message'])
        self.label.pack(pady=20)

        self.quit_button = tk.Button(master, text=self.language['quit_button'], command=self.master.quit)
        self.quit_button.pack(pady=10)