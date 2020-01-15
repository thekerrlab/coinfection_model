import pylab as pl
import sciris as sc
import shared_vars as sv


default_pars = sc.objdict({
    'popsize': 100,
    'init_prev': 0.01,
    'npts': 100,
    'beta': 0.1,
    'recov': 0.02,
    })


class HIV_Model:
    
    def __init__(self, pars, coinfection):
        self.pars = pars
        self.coinfection = coinfection
        
        self.t = 0
        self.tvec = pl.arange(pars.npts)
        self.S = pl.zeros(pars.npts)
        self.I = pl.zeros(pars.npts)
        self.R = pl.zeros(pars.npts)
        
        self.S[0] = pars.popsize*(1-pars.init_prev)
        self.I[0] = pars.popsize*pars.init_prev
    

    def step(self):
        factor = (1-sv.tb_prev) if self.coinfection else 1
        t = self.t

        dS = self.S[t]*self.pars.beta
        dR = self.I[t]*self.pars.recov*factor

        self.S[t+1] = self.S[t] - dS
        self.I[t+1] = self.I[t] + dS - dR
        self.R[t+1] = self.R[t] + dR
        self.t += 1

        sv.hiv_prev = self.I[t+1]/(self.S[t+1]+self.I[t+1]+self.R[t+1])
        print(f'HIV {t}: HIV={sv.hiv_prev:0.2f}, TB={sv.tb_prev:0.2f}')
    

    def run(self):
        for t in self.tvec[:-1]:
            self.step()
        

    def plot(self):
        pl.plot(self.tvec, self.S, label='S', lw=2)
        pl.plot(self.tvec, self.I, label='I', lw=2)
        pl.plot(self.tvec, self.R, label='R', lw=2)
        pl.title(f'"HIV" model – coinfection {self.coinfection}', fontweight='bold')
        pl.legend()
        
        
    