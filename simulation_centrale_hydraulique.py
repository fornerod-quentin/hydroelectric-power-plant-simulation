from tkinter import *
from PIL import Image, ImageTk
from time import *
import sys
import time
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.integrate import odeint
import PyPDF2
from tkinter import filedialog
import os
from itertools import chain
class main_view():
    def open_PDF(self):
        os.startfile("test.pdf")
    def SEDO(self,U,t):
        w,S = U
        J = float(self.inertia_moment)
        rho = float(self.fluid_density)
        h = float(self.dam_height)
        g = float(self.gravity)
        pole_pairs = int(self.pole_pairs)
        lg_rp = self.pw_ramp_length
        fluc = self.pw_fluctuations
        req_pow = self.required_power * 1E+6 #transform in MW
        f1 = 2 * t * t / (np.power(lg_rp,2)) #first part ramp
        f2 = 1 - 2 * np.power(((t - lg_rp) / lg_rp),2) #second part ramp
        f3 = req_pow + 2 / np.pi * np.arctan(t - lg_rp) * (fluc / 100) * req_pow * np.sin(-4 * np.pi * t)
        wc = 2 * math.pi * 50 / pole_pairs
        R = 1 / wc * math.sqrt(g * h / 2)
        v = math.sqrt(2 * g * h)
        Pdemande = req_pow
        if(t >= 0 and t <= (lg_rp / 2)):
            Pdemande = f1 * req_pow
        if(t < lg_rp and t > lg_rp / 2):
            Pdemande = f2 * req_pow
        if(t >= lg_rp):
            Pdemande = f3
        coeff_lambda = Pdemande / (wc * wc)
        #coefficients régulateur PD
        p = float(self.proportional_coeff)
        d = float(self.derivative_coeff)  
        wpoint = 1 / J * (2 * R * rho * v * S * (v - w * R) - coeff_lambda * w)
        spoint = p * (wc - w) - d * wpoint
        return[wpoint,spoint]
    def start_simulation(self):
        #parameters
        plot_resolution = 50
        t_span = np.linspace(0, self.simulation_length,self.simulation_length * plot_resolution)
        rho = float(self.fluid_density)
        h = float(self.dam_height)
        g = float(self.gravity)
        pole_pairs = int(self.pole_pairs)
        wc = 2 * math.pi * 50 / pole_pairs
        R = 1 / wc * math.sqrt(g * h / 2)
        v = math.sqrt(2 * g * h)
        #power calculation
        lg_rp = self.pw_ramp_length
        fluc = self.pw_fluctuations
        pw_req = self.required_power * 1E+6 #transform in MW
        t1 = np.linspace(0,lg_rp / 2,int(lg_rp / 2 * plot_resolution))
        t2 = np.linspace(lg_rp / 2 + 1 / plot_resolution,lg_rp,int(lg_rp / 2 * plot_resolution))
        t3 = np.linspace(lg_rp + 1 / plot_resolution,self.simulation_length,int((self.simulation_length - lg_rp) * plot_resolution))
        f1 = 2 * t1 * t1 / (np.power(lg_rp,2)) #first part ramp
        f2 = 1 - 2 * np.power(((t2 - lg_rp) / lg_rp),2) #second part ramp
        f3 = pw_req + 2 / np.pi * np.arctan(t3 - lg_rp) * (fluc / 100) * pw_req * np.sin(-4 * np.pi * t3)
        t_pw = np.array(list(chain(t1,t2,t3)))
        P1demande = f1 * pw_req
        P2demande = f2 * pw_req
        P3demande = f3
        Pdemande = np.array(list(chain(P1demande,P2demande,P3demande)))
        w_init = 0
        s_init = 0
        Uzero = [w_init, s_init]
        solution = odeint(self.SEDO, Uzero, t_span)
        torque = (2 * R * rho * solution[:,1] * v * (v - R * solution[:,0]))
        pw_turbine = torque * solution[:,0]
        coeff_lambda = Pdemande / (wc * wc)
        pw_req = coeff_lambda * solution[:,0] * solution[:,0]
        diam = np.sqrt(solution[:, 1] / (math.pi)) * 2 * 100
        fig = plt.figure()
        gs = gridspec.GridSpec(5, 1, wspace=0.25, hspace=1) # 3x1 grid
        ax0 = fig.add_subplot(gs[0, 0]) # first row
        ax0.set_title('Vitesse rotation de la turbine')
        ax0.set_xlabel('t [s]')
        ax0.set_ylabel('[tr/s] ')
        ax1 = fig.add_subplot(gs[1, 0]) # second row
        ax1.set_title('''Diamètre d'ouverture de la vanne''')
        ax1.set_xlabel('t [s]')
        ax1.set_ylabel('[cm]')
        ax2 = fig.add_subplot(gs[2, 0]) # third row
        ax2.set_title('Puissances')
        ax2.set_xlabel('t [s]')
        ax2.set_ylabel('[MW]')
        ax3 = fig.add_subplot(gs[3, 0]) # third row
        ax3.set_title('Rendement')
        ax3.set_xlabel('t [s]')
        ax3.set_ylabel('[%]')
        ax4 = fig.add_subplot(gs[4, 0]) # third row
        ax4.set_title('Couple')
        ax4.set_xlabel('t [s]')
        ax4.set_ylabel('[Nm]')
        ax0.plot(t_span, solution[:, 0] / (2 * math.pi), label='w')
        ax1.plot(t_span, diam, label='s')
        ax2.plot(t_pw,Pdemande / 1E+6,t_pw,pw_turbine / 1E+6,t_pw,pw_req / 1E+6)
        ax3.plot(t_pw,pw_req/Pdemande*100, label='rendement')
        ax4.plot(t_pw,torque)
        ax0.grid()
        ax1.grid()
        ax2.grid()
        ax3.grid()
        ax4.grid()
        ax0.legend()
        ax1.legend()
        ax2.legend(['P demande','P turbine','P fournie'])
        ax3.legend()
        ax4.legend()
        plt.show()
    def show_settings(self):
        def update_parameters():
            #simulation length
            if abs(int(sim_duration_entry.get())) > 0:
                self.simulation_length = abs(int(sim_duration_entry.get()))
                self.sim_duration_entry.config(text=self.simulation_length)
            else:
                self.simulation_length = 2 #2 seconds as minimal time
                self.sim_duration_entry.config(text=self.simulation_length)
            #dam height
            if abs(int(dam_height_entry.get()))>100:
                self.dam_height = abs(int(dam_height_entry.get()))
                self.dam_height_value_init.config(text =self.dam_height)
            else:
                self.dam_height = 100
                self.dam_height_value_init.config(text =self.dam_height)
            #power ramp duration
            if abs(float(pw_rise_time_entry.get())) > 0 and abs(float(pw_rise_time_entry.get())) < self.simulation_length:
                self.pw_ramp_length = abs(float(pw_rise_time_entry.get()))
                self.pw_rise_time_entry.config(text=self.pw_ramp_length)
            else:
                self.pw_ramp_length = 1 #1 second ramp as minimal time
                self.pw_rise_time_entry.config(text=self.pw_ramp_length)
            #power fluctuations
            if abs(int(pw_fluctuations_entry.get()))<=100:
                self.pw_fluctuations = abs(int(pw_fluctuations_entry.get()))
                self.pw_fluctuations_entry.config(text=self.pw_fluctuations)
            else:
                self.pw_fluctuations = 10
                self.pw_fluctuations_entry.config(text=self.pw_fluctuations)
            #power required
            if abs(int(pw_required_entry.get()))>0:
                self.required_power = abs(int(pw_required_entry.get()))
                self.pw_required_entry.config(text=self.required_power)
            else:
                self.required_power = 200 #minimal power [MW]
                self.pw_required_entry.config(text=self.required_power)
            #proportional coeff
            self.proportional_coeff = abs(float(prop_coeff_entry.get()))
            self.prop_coeff_entry.config(text=self.proportional_coeff)
            #derivative coeff
            self.derivative_coeff = abs(float(deriv_coeff_entry.get()))
            self.deriv_coeff_entry.config(text=self.derivative_coeff)
            #quit and destroy window
            settings_window.destroy()
            settings_window.update()
        def quit_settings():
            settings_window.destroy()
        #constants
        label_width = 20
        value_width = 10
        button_width = 40
        BG = "#FFFFFF"
        settings_window = Toplevel(self.root)
        settings_window.title('Paramètres')
        settings_frame = Frame(settings_window,bg=BG,padx=10,pady=10)
        settings_frame.grid(column=0,row=0)
        #simulation length
        sim_duration_label = Label(settings_frame,text="Durée de la simulation [s]", anchor="w",width=label_width,bg=BG)
        sim_duration_label.grid(column=0,row=0)
        sim_duration_entry = Entry(settings_frame,width=value_width,bd=1)
        sim_duration_entry.grid(column=1,row=0)
        sim_duration_entry.insert(0,self.simulation_length)
        #dam height
        dam_height_label = Label(settings_frame,text="Hauteur du barrage [m]", anchor="w",width=label_width,bg=BG)
        dam_height_label.grid(row=1,column=0)
        dam_height_entry = Entry(settings_frame,width=value_width,bd=1)
        dam_height_entry.grid(column=1,row=1)
        dam_height_entry.insert(0,self.dam_height)
        #power ramp time
        pw_rise_time_label = Label(settings_frame,text="Temps rampe puissance [s]", anchor="w",width=label_width,bg=BG)
        pw_rise_time_label.grid(row=2,column=0)
        pw_rise_time_entry = Entry(settings_frame,width=value_width,bd=1)
        pw_rise_time_entry.grid(column=1,row=2)
        pw_rise_time_entry.insert(0,self.pw_ramp_length)
        #power fluctuations
        pw_fluctuations_label = Label(settings_frame,text="Fluctuation de la puissance demandée [%]", anchor="w",width=label_width,bg=BG)
        pw_fluctuations_label.grid(row=3,column=0)
        pw_fluctuations_entry = Entry(settings_frame,width=value_width,bd=1)
        pw_fluctuations_entry.grid(column=1,row=3)
        pw_fluctuations_entry.insert(0,self.pw_fluctuations)
        #electricity required
        pw_required_label = Label(settings_frame,text="Puissance demandée [MW]", anchor="w",width=label_width,bg=BG)
        pw_required_label.grid(row=4,column=0)
        pw_required_entry = Entry(settings_frame,width=value_width,bd=1)
        pw_required_entry.grid(column=1,row=4)
        pw_required_entry.insert(0,self.required_power)# power in MW
        #proportional coeff
        prop_coeff_label = Label(settings_frame,text="Coeff proportionnel", anchor="w",width=label_width,bg=BG)
        prop_coeff_label.grid(row=5,column=0)
        prop_coeff_entry = Entry(settings_frame,width=value_width,bd=1)
        prop_coeff_entry.grid(column=1,row=5)
        prop_coeff_entry.insert(0,self.proportional_coeff)
        #derivative coeff
        deriv_coeff_label = Label(settings_frame,text="Coeff dérivateur", anchor="w",width=label_width,bg=BG)
        deriv_coeff_label.grid(row=6,column=0)
        deriv_coeff_entry = Entry(settings_frame,width=value_width,bd=1)
        deriv_coeff_entry.grid(column=1,row=6)
        deriv_coeff_entry.insert(0,self.derivative_coeff)
        #validation/interruption of the operation
        confirm_button = Button(settings_frame,text="Valider",command=lambda:update_parameters(),bg=BG,pady = 15,padx=button_width)
        confirm_button.grid(column=0,row=7)
        quit_button = Button(settings_frame,text="Annuler",command=lambda:quit_settings(),bg=BG,pady = 15,padx=button_width)
        quit_button.grid(column=1,row=7)
    def graphics_generator_window(self):
        BG = "#FFFFFF"
        graphics_window = Toplevel(self.root,bg=BG)
        graphics_window.title('Génération des graphiques')
        graphics_window.geometry("500x500")    
        graphics_frame = Frame(graphics_window,bg=BG)
        graphics_frame.grid(column=0,row=0)
        graphics_label = Label(graphics_frame,border=0, text= "Choisissez le graphique à générer",bg=BG)
        graphics_label.grid(column=0,row=0,columnspan=2)
        #buttons row 1
        power_button = Button(graphics_frame,text="Puissances",bg=BG,padx=10,pady=10)
        power_button.grid(row=1,column=0,padx=20)
        diameter_button = Button(graphics_frame,text="Diamètre vanne",bg=BG,padx=10,pady=10)
        diameter_button.grid(row=1,column=1,padx=20)
        #buttons row 2
        speed_button = Button(graphics_frame,text="Vitesse de la turbine",bg=BG,padx=10,pady=10)
        speed_button.grid(row=2,column=0,padx=20)
        volume_button = Button(graphics_frame,text="Volume du lac de rétention",bg=BG,padx=10,pady=10)
        volume_button.grid(row=2,column=1,padx=20)
        #exit button
        exit_button = Button(graphics_window,text="Quitter",command=graphics_window.destroy,bg=BG,padx=10,pady=10)
        exit_button.grid(row=3,column=0,sticky=W + E,columnspan=2)
    def show_turbine_infos(self):
        BG = "#FFFFFF"
        turbine_window = Toplevel(self.root,bg=BG)
        turbine_window.title('Turbine')
        turbine_window.geometry("800x800")
        pelton_icon = PhotoImage(file="Icons/turbine_pelton.png").subsample(1, 1)
        francis_icon = PhotoImage(file="Icons/turbine_francis.png").subsample(1, 1)
        kaplan_icon = PhotoImage(file="Icons/turbine_kaplan.png").subsample(1, 1)
        turbine_image_list = [pelton_icon,francis_icon,kaplan_icon]
        self.img_index = 0
        def scroll_left():
            if (self.img_index >= 1):
                turbine_background.configure(image=turbine_image_list[self.img_index - 1])
                self.img_index-=1                
        def scroll_right():
            if (self.img_index <= 1):
                turbine_background.configure(image=turbine_image_list[self.img_index + 1])
                self.img_index += 1
        turbine_frame = Frame(turbine_window,bg=BG)
        turbine_frame.grid(column=0,row=0)
        turbine_background = Label(turbine_window,image=pelton_icon,border=0)
        turbine_background.image = pelton_icon
        turbine_background.grid(column=0,row=1,columnspan=2)
        #left/right/quit buttons
        left_button = Button(turbine_window,text="<",command=scroll_left,font=("Helvetica", 20),padx=10,bg=BG)
        left_button.grid(row=2,column=0)
        right_button = Button(turbine_window,text=">",command=scroll_right,font=("Helvetica", 20),padx=10,bg=BG)
        right_button.grid(row=2,column=1)
        exit_button = Button(turbine_window,text="Quitter",command=turbine_window.destroy,bg=BG,padx=10,pady=10)
        exit_button.grid(row=3,column=0,sticky=W + E,columnspan=2)
    def show_barrage_infos(self):
        BG = "#FFFFFF"
        dam_window = Toplevel(self.root,bg=BG)
        dam_window.title('Barrage')
        poids_icon = PhotoImage(file="Icons/barrage_poids.png").subsample(1, 1)
        voute_icon = PhotoImage(file="Icons/barrage_voute.png").subsample(1, 1)
        contreforts_icon = PhotoImage(file="Icons/barrage_contreforts.png").subsample(1, 1)
        dam_image_list = [poids_icon,voute_icon,contreforts_icon]  
        self.img_index = 0
        def scroll_left():
            if (self.img_index >= 1):
                dam_background.configure(image=dam_image_list[self.img_index - 1])
                self.img_index-=1                
        def scroll_right():
            if (self.img_index <= 1):
                dam_background.configure(image=dam_image_list[self.img_index + 1])
                self.img_index += 1
        dam_frame = Frame(dam_window,bg=BG)
        dam_frame.grid(column=0,row=0)
        dam_background = Label(dam_window,image=poids_icon,border=0)
        dam_background.image = poids_icon
        dam_background.grid(column=0,row=1,columnspan=2)
        #left/right/quit buttons
        left_button = Button(dam_window,text="<",command=scroll_left,font=("Helvetica", 20),padx=10,bg=BG)
        left_button.grid(row=2,column=0)
        right_button = Button(dam_window,text=">",command=scroll_right,font=("Helvetica", 20),padx=10,bg=BG)
        right_button.grid(row=2,column=1)
        exit_button = Button(dam_window,text="Quitter",command=dam_window.destroy,bg=BG,padx=10,pady=10)
        exit_button.grid(row=3,column=0,sticky=W + E,columnspan=2)
    def show_vanne_infos(self):
        BG = "#FFFFFF"
        vanne_window = Toplevel(self.root,bg=BG)
        vanne_window.title('Vanne')  
        vanne_icon = PhotoImage(file="Icons/vanne.png").subsample(1, 1)
        vanne_frame = Frame(vanne_window,bg=BG)
        vanne_frame.grid(column=0,row=0)
        vanne_background = Label(vanne_window,image=vanne_icon,border=0)
        vanne_background.image = vanne_icon
        vanne_background.grid(column=0,row=1,columnspan=2)
        #quit button
        exit_button = Button(vanne_window,text="Quitter",command=vanne_window.destroy,bg=BG,padx=10,pady=10)
        exit_button.grid(row=3,column=0,sticky=W + E,columnspan=2)
    def show_alternateur_infos(self):
        BG = "#FFFFFF"
        alternateur_window = Toplevel(self.root,bg=BG)
        alternateur_window.title('Alternateur')   
        alternateur_icon = PhotoImage(file="Icons/alternateur.png").subsample(1, 1)
        alternateur_frame = Frame(alternateur_window,bg=BG)
        alternateur_frame.grid(column=0,row=0)
        alternateur_background = Label(alternateur_window,image=alternateur_icon,border=0)
        alternateur_background.image = alternateur_icon
        alternateur_background.grid(column=0,row=1,columnspan=2)
        #quit button
        exit_button = Button(alternateur_window,text="Quitter",command=alternateur_window.destroy,bg=BG,padx=10,pady=10)
        exit_button.grid(row=3,column=0,sticky=W + E,columnspan=2)
    def show_transformateur_infos(self):
        BG = "#FFFFFF"
        transformateur_window = Toplevel(self.root,bg=BG)
        transformateur_window.title('Transformateur')
        transfo_icon = PhotoImage(file="Icons/transformateur.png").subsample(1, 1)
        transformateur_frame = Frame(transformateur_window,bg=BG)
        transformateur_frame.grid(column=0,row=0)
        transformateur_background = Label(transformateur_window,image=transfo_icon,border=0)
        transformateur_background.image = transfo_icon
        transformateur_background.grid(column=0,row=1,columnspan=2)
        #quit button
        exit_button = Button(transformateur_window,text="Quitter",command=transformateur_window.destroy,bg=BG,padx=10,pady=10)
        exit_button.grid(row=3,column=0,sticky=W + E,columnspan=2)
    def __init__(self):
        #variables
        BG = "#FFFFFF"
        self.enableApplication = 1
        #creating root window in fullscreen
        self.root = Tk()
        self.root.configure(background=BG)
        self.root.title('''Simulation d'une centrale hydroélectrique''')
        #self.root.geometry("1280x800")
        self.root.attributes("-fullscreen", True)
        self.root.bind("<Escape>", lambda event: self.root.attributes("-fullscreen", False))
        #create instance image for Toplevel button (has to be there)
        self.top_frame = Frame(self.root,bg="white")
        self.top_frame.grid(row=0,column=0,sticky=W + E)
        dam_icon = PhotoImage(file="Icons/dam_picture.png").subsample(1, 1)
        self.dam_background = Label(self.top_frame,image=dam_icon,border=0)
        self.dam_background.grid(column=0,row=0)
        #initial values
        self.pole_pairs = 7
        self.dam_height = 1880
        self.valve_opening = 0
        self.inertia_moment = 100000
        self.fluid_density = 1000
        self.gravity = 9.81
        self.required_power = 423
        self.proportional_coeff = 0.05
        self.derivative_coeff = 0.04
        self.pw_ramp_length = 2 #time in seconds
        self.pw_fluctuations = 10 #power fluctuations (percentage)
        self.simulation_length = 5
        #constants
        label_width = 26
        value_width = 8
        authors_text = ("HEIG-VD, PHY2 2020 "
        "Fornerod Quentin, Hugi Samuel, Rodriguez Michaël, Vieux Adriano")
        #parameters frame
        self.side_frame = Frame(self.root,bg="white")
        self.side_frame.grid(row=0,column=1)
        #simulation duration
        self.sim_duration_label = Label(self.side_frame,text="Durée simulation [s]", anchor="w",width=label_width,bg=BG)
        self.sim_duration_label.grid(column=0,row=0)
        self.sim_duration_entry = Label(self.side_frame,text=self.simulation_length,anchor="w",width=value_width,borderwidth=2,bg=BG, relief="groove")
        self.sim_duration_entry.grid(column=1,row=0)
        #dam height
        self.dam_height_label = Label(self.side_frame,text="Hauteur du barrage [m]", anchor="w",width=label_width,bg=BG)
        self.dam_height_label.grid(row=1,column=0)
        self.dam_height_value_init = Label(self.side_frame,text=self.dam_height,anchor="w",width=value_width,bd=1,bg=BG, borderwidth=2,relief="groove")
        self.dam_height_value_init.grid(column=1,row=1)
        #power ramp duration
        self.pw_rise_time_label = Label(self.side_frame,text="Rampe puissance [s]", anchor="w",width=label_width,bg=BG)
        self.pw_rise_time_label.grid(row=2,column=0)
        self.pw_rise_time_entry = Label(self.side_frame,text=self.pw_ramp_length,anchor="w",width=value_width,bd=1,bg=BG, borderwidth=2,relief="groove")
        self.pw_rise_time_entry.grid(column=1,row=2)
        #power fluctuations
        self.pw_fluctuations_label = Label(self.side_frame,text="Fluctuations de la puissance [%]", anchor="w",width=label_width,bg=BG)
        self.pw_fluctuations_label.grid(row=3,column=0)
        self.pw_fluctuations_entry = Label(self.side_frame,text=self.pw_fluctuations,anchor="w",width=value_width,bd=1,bg=BG, borderwidth=2,relief="groove")
        self.pw_fluctuations_entry.grid(column=1,row=3)
        #power required
        self.pw_required_label = Label(self.side_frame,text="Puissance requise [MW]", anchor="w",width=label_width,bg=BG)
        self.pw_required_label.grid(row=4,column=0)
        self.pw_required_entry = Label(self.side_frame,text=self.required_power,anchor="w",width=value_width,bd=1,bg=BG, borderwidth=2,relief="groove")
        self.pw_required_entry.grid(column=1,row=4)
        #proportional coefficient
        self.prop_coeff_label = Label(self.side_frame,text="Coeff proportionnel", anchor="w",width=label_width,bg=BG)
        self.prop_coeff_label.grid(row=5,column=0)
        self.prop_coeff_entry = Label(self.side_frame,text=self.proportional_coeff,anchor="w",width=value_width,bd=1,bg=BG, borderwidth=2,relief="groove")
        self.prop_coeff_entry.grid(column=1,row=5)
        #integral coefficient
        self.deriv_coeff_label = Label(self.side_frame,text="Coeff dérivateur", anchor="w",width=label_width,bg=BG)
        self.deriv_coeff_label.grid(row=6,column=0)
        self.deriv_coeff_entry = Label(self.side_frame,text=self.derivative_coeff,anchor="w",width=value_width,bd=1,bg=BG, borderwidth=2,relief="groove")
        self.deriv_coeff_entry.grid(column=1,row=6)
        #open settings window
        self.settings_button = Button(self.side_frame,text="Paramètres",command=lambda:self.show_settings(),pady=15,bg=BG)
        self.settings_button.grid(column=0,row=7,sticky=W + E)
        #bottom frame
        self.bottom_frame = Frame(self.root,bg=BG)
        self.bottom_frame.grid(row=1,column=0)
        #start graphics calculation
        self.start_button = Button(self.bottom_frame,text="Générer les graphiques",command=lambda:self.start_simulation(),bg=BG,pady = 15)
        self.start_button.grid(column=0,row=0)
        #open turbine window
        self.turbine_button = Button(self.bottom_frame,text="Infos turbine",command=lambda:self.show_turbine_infos(),pady=15,bg=BG)
        self.turbine_button.grid(column=1,row=0)
        #open barrage window
        self.barrage_button = Button(self.bottom_frame,text="Infos barrage",command=lambda:self.show_barrage_infos(),bg=BG,pady = 15)
        self.barrage_button.grid(column=2,row=0)
        #open vanne_régulation window
        self.vanne_button = Button(self.bottom_frame,text="Infos vanne",command=lambda:self.show_vanne_infos(),bg=BG,pady = 15)
        self.vanne_button.grid(column=3,row=0)
        #open alternateur window
        self.alternateur_button = Button(self.bottom_frame,text="Infos alternateur",command=lambda:self.show_alternateur_infos(),bg=BG,pady = 15)
        self.alternateur_button.grid(column=4,row=0)
        #open transformateur window
        self.transformateur_button = Button(self.bottom_frame,text="Infos transformateur",command=lambda:self.show_transformateur_infos(),bg=BG,pady = 15)
        self.transformateur_button.grid(column=5,row=0)
        #PDF rapport
        self.pdf_button = Button(self.bottom_frame,text="Rapport en PDF",command=lambda:self.open_PDF(),pady=15,bg=BG)
        self.pdf_button.grid(column=7,row=0)
        #exit button
        self.exit_button = Button(self.bottom_frame,text="Quitter",command=self.quit,pady=15,bg=BG)
        self.exit_button.grid(column=8,row=0)
        #copyright/authors
        end_frame = Frame(self.root,bg=BG)
        end_frame.grid(row=2,column=0)
        self.authors_label = Label(end_frame,text=authors_text,font=("Helvetica", 8), bg=BG)
        self.authors_label.grid(column=0,row=1,pady=10)
        #launch while loop
        self.mainProgram()
    def toggle_start(self):
        if(self.start_state):
            self.start_state = not(self.start_state)            
    def mainProgram(self):
        while True and self.enableApplication:
            try:
                self.root.update()
                self.root.update_idletasks()
            except:
                break
        return
    def quit(self):
       self.enableApplication = 0
       sleep(1)
       self.root.destroy()
       sys.exit()
def main():
    view1 = main_view()
if __name__ == '__main__':
    main()