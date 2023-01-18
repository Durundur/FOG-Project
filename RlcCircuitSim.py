import numpy as np
import matplotlib.pyplot as plt
import PySpice.Logging.Logging as Logging
from PySpice.Spice.Netlist import Circuit
from PySpice.Unit import *
from matplotlib.widgets import Slider, RadioButtons
from sympy.physics.quantum.circuitplot import matplotlib


class RlcCircuitSim:

    def __init__(self, resistance, inductance, capacitance, omega, source_signal, amplitude):
        self.resistance = resistance
        self.inductance = inductance
        self.capacitance = capacitance
        self.omega = omega
        self.source_signal = source_signal
        self.amplitude = amplitude
        self.freq = self.omega / (2 * np.pi)
        self.create_gui()
        self.update_graphs()
        matplotlib.use('Qt5Agg')
        plt.show()

    def handle_sliders_change(self, Null):
        self.resistance = self.r_slider.val
        self.inductance = self.l_slider.val
        self.capacitance = self.c_slider.val
        self.amplitude = self.a_slider.val
        self.omega = self.omega_slider.val
        self.freq = self.omega / (2 * np.pi)
        self.update_graphs()

    def update_graphs(self):
        self.r_graph2.cla()
        self.l_graph2.cla()
        self.c_graph2.cla()
        self.r_graph.cla()
        self.l_graph.cla()
        self.c_graph.cla()
        self.source_graph.cla()
        self.voltages_graph.cla()
        self.r_graph.grid()
        self.l_graph.grid()
        self.c_graph.grid()
        self.source_graph.grid()
        self.voltages_graph.grid()
        analysis = self.simulate_circuit(self.amplitude, self.freq, self.resistance, self.capacitance, self.inductance,self.source_signal)
        time = np.array(analysis.time)
        source_voltage = np.array(analysis["n1"])
        r_voltage = np.zeros_like(time)
        l_voltage = np.zeros_like(time)
        c_voltage = np.zeros_like(time)
        if self.resistance != 0:
            if self.capacitance == 0 and self.inductance == 0:
                r_voltage = np.array(analysis["n1"])
            if self.inductance != 0:
                if self.capacitance != 0:
                    r_voltage = np.array(analysis["n1"]) - np.array(analysis["n2"])
                    l_voltage = np.array(analysis["n2"]) - np.array(analysis["n3"])
                    c_voltage = np.array(analysis["n3"])
                else:
                    r_voltage = np.array(analysis["n1"]) - np.array(analysis["n2"])
                    l_voltage = np.array(analysis["n2"])
            else:
                if self.capacitance != 0:
                    c_voltage = np.array(analysis["n2"])
                else:
                    r_voltage = np.array(analysis["n1"])
        else:
            if self.inductance != 0:
                if self.capacitance != 0:
                    l_voltage = np.array(analysis["n1"]) - np.array(analysis["n2"])
                    c_voltage = np.array(analysis["n2"])
                else:
                    l_voltage = np.array(analysis["n1"])
            else:
                if self.capacitance != 0:
                    c_voltage = np.array(analysis["n1"])

        r_current = np.zeros_like(time)
        l_current = np.zeros_like(time)
        c_current = np.zeros_like(time)
        if self.resistance != 0:
            r_current = np.array(analysis["vr1_minus"])
        if self.inductance != 0:
            l_current = np.array(analysis["vl1_minus"])
        if self.capacitance != 0:
            c_current = np.array(analysis["vc1_minus"])
        self.source_graph.set_title("Sygnał źródła")
        self.source_graph.set_xlabel('t [s]')
        self.source_graph.set_ylabel('U [V]', labelpad=-1)
        self.source_graph.plot(time, source_voltage, color='blue')
        self.r_graph.set_title('Rezystor')
        self.r_graph.set_xlabel('t [s]')
        self.r_graph.set_ylabel('U [V]', labelpad=-1, color='orange')
        self.r_graph.plot(time, r_voltage, color='orange')
        self.r_graph2.set_ylabel('I [A]', labelpad=-1, color='gray')
        self.r_graph2.plot(time, r_current, color='gray')
        self.l_graph.set_title('Cewka')
        self.l_graph.set_xlabel('t [s]')
        self.l_graph.set_ylabel('U [V]', labelpad=-1, color='green')
        self.l_graph.plot(time, l_voltage, color='green')
        self.l_graph2.plot(time, l_current, color='gray')
        self.l_graph2.set_ylabel('I [A]', labelpad=-1, color='gray')
        self.c_graph.set_title('Kondensator')
        self.c_graph.set_xlabel('t [s]')
        self.c_graph.set_ylabel('U [V]', labelpad=-1, color='red')
        self.c_graph.plot(time, c_voltage, color='red')
        self.c_graph2.plot(time, c_current, color='grey')
        self.c_graph2.set_ylabel('I [A]', labelpad=-1, color='gray')
        self.voltages_graph.set_xlabel('t [s]')
        self.voltages_graph.set_ylabel('U [V]', labelpad=-1)
        self.voltages_graph.plot(time, source_voltage, label=r'$U_{Źródła}$', color='blue')
        self.voltages_graph.plot(time, r_voltage, label=r'$U_{Rezystora}$', color='orange')
        self.voltages_graph.plot(time, l_voltage, label=r'$U_{Cewki}$', color='green')
        self.voltages_graph.plot(time, c_voltage, label=r'$U_{Kondensatora}$', color='red')
        self.voltages_graph.set_title("Wykresy napięć")
        self.voltages_graph.legend()
        plt.draw()

    def simulate_circuit(self, amplitude, freq, resistance, capacitance, inductance, source_signal):
        logger = Logging.setup_logging()
        circuit = Circuit('RLC')
        if source_signal == "prostokąt":
            Vac = circuit.PulseVoltageSource('input', 'n1', circuit.gnd, initial_value=-amplitude @ u_V,
                                             pulsed_value=amplitude @ u_V, pulse_width=1 / (freq * 2), period=1 / freq)
        elif source_signal == "sinus":
            Vac = circuit.SinusoidalVoltageSource('input', 'n1', circuit.gnd, amplitude=amplitude @ u_V,
                                                  frequency=freq @ u_Hz)
        elif source_signal == "cosinus":
            Vac = circuit.SinusoidalVoltageSource('input', 'n1', circuit.gnd, amplitude=amplitude @ u_V,
                                                  frequency=freq @ u_Hz, delay=-1 / (4 * freq))
        elif source_signal == "piła":
            Vac = circuit.PulseVoltageSource('input', 'n1', circuit.gnd, initial_value=-amplitude @ u_V,
                                             pulsed_value=amplitude @ u_V, fall_time=0, rise_time=1 / freq,
                                             pulse_width=1 / freq, period=1 / freq, delay_time=-1 / freq / 2)
        inserted_elements_counter = 0
        elements_array = [resistance, inductance, capacitance]
        for i in range(3):
            start_node = 'n' + str(inserted_elements_counter + 1)
            end_node = 'n' + str(inserted_elements_counter + 2)
            if i == 2 or (elements_array[i + 1] == 0 and i == 1) or (
                    i == 0 and elements_array[i + 1] == 0 and elements_array[i + 2] == 0):
                end_node = circuit.gnd
            if elements_array[i] == 0:
                continue
            if i == 0:
                circuit.R(1, start_node, end_node, resistance @ u_Ohm)
                inserted_elements_counter = inserted_elements_counter + 1
                circuit.R1.minus.add_current_probe(circuit)

            if i == 1:
                circuit.L(1, start_node, end_node, inductance @ u_mH)
                inserted_elements_counter = inserted_elements_counter + 1
                circuit.L1.minus.add_current_probe(circuit)
            if i == 2:
                circuit.C(1, start_node, end_node, capacitance @ u_uF)
                inserted_elements_counter = inserted_elements_counter + 1
                circuit.C1.minus.add_current_probe(circuit)

        circuit.Vinput.minus.add_current_probe(circuit)
        simulator = circuit.simulator()
        analysis = simulator.transient(step_time=0.0001, end_time=Vac.period * 100 + 0.2, start_time=Vac.period * 100)
        # analysis = simulator.transient(step_time=0.0001, end_time=0.1, start_time=0)
        return analysis

    def handle_source_signal_change(self, label):
        self.source_signal = label
        self.update_graphs()

    def create_gui(self):
        self.fig = plt.figure()
        self.r_graph = plt.axes([.05, .57, .25, .40])
        self.r_graph.grid()
        self.r_graph2 = self.r_graph.twinx()
        self.l_graph = plt.axes([.37, .57, .25, .40])
        self.l_graph2 = self.l_graph.twinx()
        self.c_graph = plt.axes([.69, .57, .25, .40])
        self.c_graph2 = self.c_graph.twinx()
        self.source_graph = plt.axes([.68, .08, .28, .40])
        self.voltages_graph = plt.axes([.35, .08, .28, .40])
        self.omega_axe = plt.axes([.043, .28, .25, .016])
        self.omega_slider = Slider(self.omega_axe, r'$\omega$', 1, 250, valstep=1)
        self.omega_slider.set_val(self.omega)
        self.a_axe = plt.axes([.043, .23, .25, .016])
        self.a_slider = Slider(self.a_axe, 'A [V]', 1, 150, valstep=1)
        self.a_slider.set_val(self.amplitude)
        self.c_axe = plt.axes([.043, .18, .25, .016])
        self.c_slider = Slider(self.c_axe, 'C [\u03BCF]', 0, 500, valstep=1)
        self.c_slider.set_val(self.capacitance)
        self.l_axe = plt.axes([.043, .13, .25, .016])
        self.l_slider = Slider(self.l_axe, 'L [mH]', 0, 1000, valstep=1)
        self.l_slider.set_val(self.inductance)
        self.r_axe = plt.axes([.043, .08, .25, .016])
        self.r_slider = Slider(self.r_axe, 'R [Ohm]', 0, 500, valstep=1)
        self.r_slider.set_val(self.resistance)
        self.source_signal_change_axe = plt.axes([.03, .33, .08, .15])
        self.source_signal_change_axe.set_title("Sygnał na źródle")
        self.source_signal_change_radio = RadioButtons(self.source_signal_change_axe,
                                                       ('prostokąt', 'sinus', 'cosinus', 'piła'))
        self.source_signal_change_radio.on_clicked(self.handle_source_signal_change)
        sliders = [self.r_slider, self.l_slider, self.c_slider, self.a_slider, self.omega_slider]
        for slide in sliders:
            slide.on_changed(self.handle_sliders_change)


RlcCircuitSim(0, 700, 370, 100, "prostokąt", 150)

