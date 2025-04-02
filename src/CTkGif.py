from PIL import Image
import customtkinter as ctk
import threading
import os


class CTkGif(ctk.CTkLabel):
    """ Affichage d'un gif """

    def __init__(self, master: any, path, loop=False, acceleration=1, repeat=1, size=None, stop_command=None, **kwargs):
        super().__init__(master, **kwargs)
        if acceleration <= 0:
            raise ValueError('Acceleration must be strictly positive')
        self.master = master  # master
        self.repeat = repeat  # Nombre repeat max
        self.path = path  # Chemin gif
        self.count = 0  # Nombre d'animations réalisées
        self.loop = loop  # tourne on a l'infini ?
        self.acceleration = acceleration  # facteur de ralentissement
        self.index = 0  # Frame affichée
        self.is_playing = False  # état de l'affichage
        try:
            self.gif = Image.open(path)  # Image ouverte
        except FileNotFoundError:
            print(os.listdir('data/ressources'))
            exit(1)
        self.n_frame = self.gif.n_frames  # Nombres de frame de l'animation
        self.frame_duration = self.gif.info['duration'] * 1/self.acceleration  # temps d'une frame
        self.force_stop = False
        self.stop_command = stop_command
        if size is not None:
            self.size = size
        else:
            self.size = self.gif.size
        self.configure(text='', image=ctk.CTkImage(self.gif, size=self.size))

    def update(self):  # Update l'affichage du gid
        if self.index < self.n_frame:  # Si on est pas au bout de l'animation
            if not self.force_stop:  # Si on est pas forcé de s'arrêter
                self.gif.seek(self.index)  # Frame suivante
                self.configure(image=ctk.CTkImage(self.gif, size=self.size))  # Affichage
                self.index += 1  # Indexation
                self.after(int(self.frame_duration), self.update)  # Programmation prochaine frame
            else:
                self.force_stop = False  # On passe en off
                if self.stop_command is not None:
                    self.stop_command()
        else:  # Si on est au bout
            self.index = 0  # On revient au début
            self.count += 1  # On incrémente le compteur
            if self.is_playing and (self.count < self.repeat or self.loop):  # pas d'arrêt et on recommence
                self.after(int(self.frame_duration), self.update)  # Programmation prochaine frame
            else:
                self.is_playing = False  # On passe en off
                if self.stop_command is not None:
                    self.stop_command()

    def start(self):
        """ débute l'animation si arrêtée"""
        if not self.is_playing:
            self.count = 0
            self.is_playing = True
            self.after(int(self.frame_duration), self.update)

    def stop(self, forced=False):
        """arrête l'animation brusquement si forcé, à sa fin sinon"""
        if self.is_playing:  # Si on joue
            self.is_playing = False
            self.force_stop = forced

    def toggle(self, forced=False):
        """ change le status de la lecture"""
        if self.is_playing:
            self.stop(forced=forced)
        else:
            self.start()

