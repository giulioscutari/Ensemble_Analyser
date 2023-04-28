import numpy as np
import os
from scipy.constants import c, h, electron_volt, R
import matplotlib.pyplot as plt


from ensemble_analyser.regex_parsing import regex_parsing

FACTOR_EV_NM = h*c/(10**-9*electron_volt)


class Graph:

    def __init__(self, confs, protocol, log, T, uv_ref = '', ecd_ref = '', final_lambda = 800., definition=4):
        """
        
        confs | list : whole ensemble list
        protocol | Protocol : protocol that calculates the electronic spectra 
        log : logger instance
        T | float : temperature [K]
        uv_ref | str : UV experimental graph
        ecd_ref | str : ECD experimental grpah
        final_lambda | float : last wavelength to convolute the spectra
        definition | int : number of point for the wavelength interval
        """

        self.confs = [i for i in confs if i.active]
        self.protocol = protocol
        self.log = log
        self.pop = self.calc_pop(T)
        self.log.debug(self.pop)

        self.x = np.linspace(FACTOR_EV_NM/1, FACTOR_EV_NM/(final_lambda), 10**definition) # eV x axis

        self.spectra = []
        self.filter_outputs()


        self.uv_impulses = [self.get_uv(i) for i in self.spectra]
        self.ecd_impulses = [self.get_ecd(i) for i in self.spectra]


        uv = self.calc_graph(impulses=self.uv_impulses, sigma=1/3, fname=        f"uv_protocol_{self.protocol.number}.dat", save=True)
        ecd = self.calc_graph(impulses=self.ecd_impulses, sigma=1/3, fname=        f"ecd_protocol_{self.protocol.number}.dat", save=True)

        Graph.damp_graph(f'uv_protocol_{self.protocol.number}.dat', self.x, uv)
        Graph.damp_graph(f'ecd_protocol_{self.protocol.number}.dat', self.x, ecd)
        

    def filter_outputs(self) -> None:
        """
        Read once all conformers output, keeping only graph information
        
        return None
        """

        st , en = regex_parsing[self.protocol.calculator]['start_spec'], regex_parsing[self.protocol.calculator]['end_spec']

        for i in self.confs:
            with open(os.path.join(os.getcwd(), i.folder, f'protocol_{self.protocol.number}.out')) as f:
                fl = f.read()

            sp = fl.split(st)[-1].split(en)[0]
            self.spectra.append(sp)

        return None

    def get_uv(self, spectra):
        """
        Get the impulses for the UV spectra calculation

        spectra | str : parse output file

        return | tuple(float, float) : energy and impulse tuple
        """
        graph = spectra.split(regex_parsing[self.protocol.calculator]['s_UV'])[-1].split(regex_parsing[self.protocol.calculator]['break_spectra'])[0]

        return [(FACTOR_EV_NM/float(i.strip().split()[regex_parsing[self.protocol.calculator]['idx_en_UV']]), float(i.strip().split()[regex_parsing[self.protocol.calculator]['idx_imp_UV']])) for i in graph.splitlines() if i]
    
    def get_ecd(self, spectra):
        """
        Get the impulses for the ECD spectra calculation

        spectra | str : parse output file

        return | tuple(float, float) : energy and impulse tuple
        """
        graph = spectra.split(regex_parsing[self.protocol.calculator]['s_ECD'])[-1].split(regex_parsing[self.protocol.calculator]['break_spectra'])[0]

        return [(FACTOR_EV_NM/float(i.strip().split()[regex_parsing[self.protocol.calculator]['idx_en_ECD']]), float(i.strip().split()[regex_parsing[self.protocol.calculator]['idx_imp_ECD']])) for i in graph.splitlines() if i]


    def gaussian(self, x, ev, I, sigma) -> np.array:
        """
        Create a gaussian convolution for each impulse 

        ev | float : energy of the impulse
        I | float : intensity of the impulse (Fosc for UV, R(vel) for ECD)
        sigma | float : sigma of the gaussian distribution

        I/(σ*sqrt(2π)) * exp(-1/2*((x-ev)/σ)^2)

        return | np.array : Gaussian of the impulse
        """
        return I/(sigma*np.sqrt(2*np.pi)) * np.exp(-0.5*((x-ev)/sigma)**2)

    def calc_pop(self, T):
        """
        Boltzmann population, if necessary with correction

        T | float : temperature [K]

        return | np.array (1D) : population distribution
        """
        
        # if not self.protocol.take_from:
        n , n_1 = self.protocol.number, str(int(self.protocol.number)-1)
        # else: 
        #     n , n_1 = self.protocol.number, str(self.protocol.take_from)

        # CONF do have frequency calculation before
        if self.protocol.freq:
            ens = np.array([i.get_energy for i in self.confs])
        elif self.confs[0].energies[n_1]['G']:
            # So energy is corrected: if functionals are the same, nothing change; else energy of the new function is corrected with lower frequency correction
            ens = np.array([
                i.energies[n]['E'] + (i.energies[n_1]['G'] - i.energies[n_1]['E'])
                for i in self.confs
                ])
        else: 
            raise IOError('No frequency calculation')

        ens_rel = ens - min(ens)
        bolz = np.exp((-ens_rel*4186)/(R*T))
        pop = (bolz/np.sum(bolz))
        for idx, i in enumerate(list(ens)):
            self.confs[idx]._last_energy['G'] = ens[idx]
            self.confs[idx]._last_energy['Erel'] = i
            self.confs[idx]._last_energy['Pop'] = pop[idx] * 100

        return pop


    def calc_graph(self, impulses, sigma, fname = '', save=False):
        """
        Build the UV spectra 

        impulses | list(tuple) : list of the (eV, I)s for each root of each conformer
        sigma | float : dispersion for the gaussian convolution [eV]
        fname | str : the name of the file to store the graph
        save | bool : save the .dat file for each conformer

        return | np.array (1D array) : Convoluted and weighted spectra
        """
        
        x = self.x.copy()
        y = np.zeros(x.shape)

        for idx in range(len(self.confs)):
            y_ = np.zeros(x.shape)
            for ev, I in impulses[idx]:
                y_ += self.gaussian(x, ev, I, sigma)

            if save: Graph.damp_graph(
                        fname = os.path.join(os.getcwd(), self.confs[idx].folder, fname), 
                        x = x, y = y_
                    )

            y += self.pop[idx] * y_

        return y


    @staticmethod
    def damp_graph(fname, x, y):
        """
        Damp an xy graph into file 

        fname | str : filename to store the graph
        x | np.array (1D) : x axis
        y | np.array (1D) : y axis

        return None
        """
        data = np.array([(xi, yi) for xi, yi in zip(x, y)])
        np.savetxt(fname, data)
        return 
    
    @staticmethod
    def load_graph(fname, is_ev=False):
        arr = np.loadtxt(fname)
        if not is_ev:
            return FACTOR_EV_NM/arr[:, 0], arr[:, 1]
        return arr[:, 0], arr[:, 1]