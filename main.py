import matplotlib.pyplot as plt
import numpy as np

class TurboJetSimulator:

    def __init__(self):

        pass

    def complexlogic(self,Ta=None,pa=None,TTa=None,TT3=None,TT4=None, TT5=None,
                     T6=None,pTa=None,pT2=None,pT3=None,pT4=None,pT5=None,p6=None,pcr=None,
                     ua=None,u6=None,M=None,pc=None,pb=None,nc=None,nd=None,npt=None,npc=None,
                     nt=None,massflow=None,f=None,matching=None, Thrust=None, TT4TT2=None, EPR=None,A=None, ro6=None, TSFC=None, pd=None ):

        p = {
            "Ta" : Ta,
            "pa" : pa,
            "TTa" : TTa,
            "TT3" : TT3,
            "TT4" : TT4,
            "TT5" : TT5,
            "pTa" : pTa,
            "pT2" : pT2,
            "pT3" : pT3,
            "pT4" : pT4,
            "pT5" : pT5,
            "pcr" : pcr,
            "p6" : p6,
            "T6" : T6,
            "ua" : ua,
            "u6" : u6,
            "ro6": ro6,
            "M" : M,
            "pd" : pd,
            "pc" : pc,
            "pb" : pb,
            "nc" : nc,
            "nd" : nd,
            "npc" : npc,
            "npt" : npt,
            "nt" : nt,
            "y1" : 1.4,
            "y2" : 1.33,
            "nm" : 0.99,
            "cpa" : 1005,
            "cpg" : 1147,
            "cpm" : 1076,
            "Q" : 45000000,
            "massflow" : massflow,
            "matching" : matching,
            "f" : f, 
            "A": A,
            "Thrust" : Thrust,
            "TT4TT2": TT4TT2,
            "EPR": EPR,
            "TSFC": TSFC
           
            
        }

        formulas = {
            "Ta" : lambda p: (p['TTa']/((p['pTa']/p['pa'])**((p['y1']-1)/p['y1']))) if Ta is None else None,
            "TTa": lambda p: (p['Ta']*(1 + ((p['y1']-1)/2)*p['M']**2)) if TTa is None else None,
            "pTa": lambda p: (p['pa']*(p['TTa']/p['Ta'])**(p['y1']/(p['y1']-1))) if pTa is  None else None,
            "pT2": lambda p: ((p['pa']*(1 + p['nd']*((p['y1']-1)/2)*p['M']**2)**(p['y1']/(p['y1'] - 1))) if pT2 is None and M is not None and nd is not None
                              else (p['pTa']*p['pd'] if pd is not None else None)),
            "pT3": lambda p: p['pc']*p['pT2'],
            "pT4": lambda p: p['pb']*p['pT3'],
            "TT3": lambda p: (p['TTa']*( 1 + (1/p['nc'])*(p['pc']**((p['y1']-1)/p['y1']) - 1)) if nc is not None 
                              else (p['TTa']*((p['pT3']/p['pT2'])**((p['y1']-1)/(p['npc']*p['y1']))) if npc is not None else None)),
            "TT4": lambda p: p['TT5'] + p['cpa']*(p['TT3'] - p['TTa'])/(p['cpg']*p['nm']) if TT4 is None else None,
            "TT5": lambda p: p['TT4'] - p['cpa']*(p['TT3'] - p['TTa'])/(p['cpg']*p['nm']) if TT5 is None else None,
            "pT5": lambda p: (p['pT4']*(p['TT5']/p['TT4'])**(p['y2']/((p['y2']-1)*p['npt'])) if npt is not None 
                              else (p['pT2']*p['EPR'] if EPR is not None else None)),
            "pcr": lambda p: p['pT5']*0.54,
            "T6": lambda p: p['TT5']/((p['pT5']/p['p6'])**((p['y2'] - 1)/p['y2'])),
            "p6": lambda p: (p['pa'] if p['matching']=='yes' else (p['pcr'] if p['pcr'] > p['pa'] else p['pa'])),
            "u6": lambda p: (2*p['cpg']*(p['TT5'] - p['T6']))**0.5,
            "ro6": lambda p: 1000*p['p6']/(287*p['T6']),
            "massflow": lambda p: p['ro6']*p['u6']*p['A'] if massflow is None else None,
            "A": lambda p: p['massflow']/(p['ro6']*p['u6']) if A is None else None,
            "ua": lambda p: (2*p['cpa']*(p['TTa'] - p['Ta']))**0.5,
            "f": lambda p: p['cpm']*(p['TT4'] - p['TT3'])/(p['Q']*p['nt']),
            "Thrust": lambda p: p['massflow']*((1 + p['f'])*p['u6'] - p['ua']) + p['A']*(p['p6']*1000 - p['pa']*1000),
            "TT4TT2": lambda p: p['TT4']/p['TTa'],
            "EPR": lambda p: p['pT5']/p['pT2'],
            "TSFC": lambda p: p['f']*p['massflow']*3600/p['Thrust']
        }

        attempts = 0

        while any(p[key] is None for key in p) and attempts < 5:

            

            

            for key in p:

                

                if key in formulas:

                    print(f"\n trying to compute key: {key}")

                    try:

                        result = formulas[key](p)

    

                        if result is not None:

                            try:

                               p[key] = result
                               print(p[key])

                            except Exception as e:

                                print(f"Error for {key}:{e}")

                                pass

                    except Exception as e:

                        print(f"Error for {key}:{e}")
                        pass
              
            attempts = attempts + 1

        TT3corrected = p['TTa']*((p['pT3']/p['pT2'])**((p['y1']-1)/p['y1']))

        TTdx = np.full(100,p['TTa'])
        pTdx = np.linspace(p['pTa'], p['pT2'], 100)

        TTcx = np.linspace(p['TTa'],TT3corrected, 100)
        pTcx = p['pT2'] *(TTcx / p['TTa']) ** (p['y1']/(p['y1'] - 1))

        TTtx = np.linspace(p['TT4'], p['TT5'], 100)
        pTtx = p['pT4'] * (TTtx / p['TT4'] ) ** (p['y2']/(p['y2'] - 1 ))

        TTbx = np.linspace(p['TT3'], p['TT4'], 100)
        pTbx = np.linspace(p['pT3'], p['pT4'], 100)

        TTnx = np.full(100,p['TT5'])
        pTnx = np.linspace(p['pT5'], p['p6'],100)

        


        X = np.concatenate([TTdx, TTcx, TTbx, TTtx,TTnx])
        P = np.concatenate([pTdx, pTcx, pTbx, pTtx, pTnx])

        plt.plot(X,P)
        plt.title("Total Temperature across the turbojet")
        plt.xlabel("Temperature K")
        plt.ylabel("Pressure kPa")
        plt.grid(True)

        labels = ['TTa', 'TT3', 'TT4', 'TT5']

        positions = [ 100, 200, 300 ,400]

        for label , i in zip(labels, positions):
            plt.annotate(label, (X[i], P[i]),
                         textcoords = "offset points", xytext=(0,-15),
                         ha='center', fontsize=11, color='red')

    

        return p
            
            
            
