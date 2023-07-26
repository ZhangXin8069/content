import numpy as np
import gvar as gv
import data_tools as aly
import matplotlib.pyplot as plt
import lsqfit


def main():
    gagv_aly = GVGA_ALY()
    #gagv_aly.read_3pt(curr=0)#read 3pt data
    #gagv_aly.read_2pt()#read 2pt data
    gagv_aly.ratio_32pt(curr='gv')#fit 3pt/2pt ratio
    #gagv_aly.ratio_gav()#ratio ga/gv

    

class GVGA_ALY(object):
    def read_3pt(self, curr=0):
        data_3pt = {}
        for tseq in [3,4,5,6,7]:
            data = np.load("./data/data_3pt_tseq"+str(tseq)+".npy")
            data = data[:,:,0:(tseq+1),-2:]
            #print(data.shape)#(conf, ga/gv, t, re/im)
            data = aly.jackknife_D4(data)#resample
            data_gv = data[:,curr,:,0]
            #print(data_gv.shape)#(conf, t)

            combine_gv = np.sqrt(data_gv[:,:] * data_gv[:,::-1])#Symmetrical operations for 3pt
            gva_cv = np.mean(combine_gv, axis=0)
            gva_err = np.std(combine_gv, axis=0)*np.sqrt(50)
            data_3pt['tseq='+str(tseq)] = gv.gvar(gva_cv, gva_err)
            print(data_3pt['tseq='+str(tseq)].shape)
        return(data_3pt)#return key{tseq=3,4,5...}, each have number shape (4,5,6...)

    def read_2pt(self):
        data = np.load("./data/data_2pt.npy")
        data = data[:,0,:,-2:]
        #print(data.shape)#(conf, t, re/im)

        data = aly.jackknife_D3(data)
        data = data[:,:,0]#shape(conf,t)
        ha_cv = np.mean(data, axis=0)
        ha_err = np.std(data, axis=0)*np.sqrt(50)
        print(gv.gvar(ha_cv, ha_err).shape)
        return(gv.gvar(ha_cv, ha_err))#return shape(t) 
     
    def ratio_32pt(self, curr, fig=True):
        data_3pt ={}
        if (curr == "gv"): 
            data_3pt = self.read_3pt(curr=0)
        elif (curr == "ga"): 
            data_3pt = self.read_3pt(curr=1)
        data_2pt = self.read_2pt()

        count = 0
        ratio_32 ={}
        for tseq in [3,4,5,6]:###plot data
            data_2ptseq = data_2pt[tseq]
            ratio_32['tseq='+str(tseq)] = data_3pt['tseq='+str(tseq)]/data_2ptseq

            cv = [val.mean for val in ratio_32['tseq='+str(tseq)]]
            err = [val.sdev for val in ratio_32['tseq='+str(tseq)]]
            xt = np.arange(0,tseq+1)-tseq/2
            if fig==True:   plt.errorbar(xt, cv, err, capsize=4.0,fmt='D', ms=3.0,color=col[count], label = '3pt/2pt, tseq='+str(tseq))
            count +=1

        def fit_function(x, paras):
            ans = {}# ratio 
            for tseq in [3,4,5,6]:
                ans['tseq='+str(tseq)] =  paras['p3_C'] * (1 + paras['p3_C1'] * (np.exp(-paras['DeltaE']*x['tseq='+str(tseq)]) \
                                                                        + np.exp(-paras['DeltaE']*(tseq - x['tseq='+str(tseq)])) ))
            ans['2pt'] = paras['p2_C'] * (1 + paras['p2_C1'] * np.exp(-paras['DeltaE'] * x['t']))
            return ans
        
        x_fit = {}; y_fit = {}
        x_fit['t'] = np.arange(2,8)
        y_fit['2pt'] = data_2pt[2:8]
        for tseq in [3,4,5,6]:
            x_fit['tseq='+str(tseq)] = np.arange(0, tseq+1)
            y_fit['tseq='+str(tseq)] = ratio_32['tseq='+str(tseq)]

        prior = gv.gvar(dict(p2_C='10000.0(100000)', p2_C1='1.0(100000)', p3_C='1.0(10)', p3_C1='1.0(100)', DeltaE='1.0(10)'))
        fit = lsqfit.nonlinear_fit(data=(x_fit, y_fit), svdcut=1e-8, prior=prior, fcn=fit_function)
        print(fit.format(maxline=True))     

        re_data = fit_function(x_fit, fit.p)#reproduce fitted data
        print(re_data)
        xx = np.arange(0,5+1)-5/2
        if fig==True:   plt.fill_between(xx, fit.p['p3_C'].mean-fit.p['p3_C'].sdev, fit.p['p3_C'].mean+fit.p['p3_C'].sdev, color='red', alpha=0.3)
        
        count = 0
        for tseq in [3,4,5,6]:#plot fitted data
            cv = np.array([val.mean for val in re_data['tseq='+str(tseq)]])
            err = np.array([val.sdev for val in re_data['tseq='+str(tseq)]])
            xx = np.arange(0,tseq+1)-tseq/2
            if fig==True:   plt.fill_between(xx, cv+err, cv-err, color=col[count], alpha=0.3)
            count +=1

        if fig==True:   
            #plt.axhline(0,linestyle='--',color='black') 
            plt.xlabel('t-tseq/2')
            if curr == "gv":    plt.ylim(0.6,0.9)
            if curr == "ga":    plt.ylim(0.8,1.15)
            plt.title('3pt/2pt_'+curr)
            plt.legend();   
            plt.savefig('./3pt_2pt_'+curr+'.png')
            plt.show()
        return  fit.p['p3_C']

    def ratio_gav(self):
        data_gv = self.read_3pt(curr=0)
        data_ga = self.read_3pt(curr=1)
        print(data_gv.keys())

        gv = self.ratio_32pt(curr='gv', fig=False)#load fittet gv
        ga = self.ratio_32pt(curr='ga', fig=False)#load fittet ga
        ga_r_gv = ga/gv
        print('Final ga='+str(ga_r_gv))
        count = 0
        for tseq in [3,4,5,6]:#plot ga/gv for each 3pt/2pt
            ratio = data_ga['tseq='+str(tseq)]/data_gv['tseq='+str(tseq)]
            print(ratio)
            cv = [val.mean for val in ratio]
            err = [val.sdev for val in ratio]
            xt = np.arange(0,tseq+1)-tseq/2
            plt.errorbar(xt, cv, err, capsize=4.0,fmt='D', ms=3.0,color=col[count], label = 'ga/gv, tseq='+str(tseq))
            count +=1
        #plot fitted_ga/fitted_gv
        plt.fill_between(np.arange(-3,4), ga_r_gv.mean-ga_r_gv.sdev, ga_r_gv.mean+ga_r_gv.sdev, color='red', alpha=0.3)
        #plt.axhline(0,linestyle='--',color='black') 
        plt.xlabel('t-tseq/2')
        plt.title('ga/gv')
        plt.ylim(1,1.5)
        plt.legend();   
        plt.savefig('./ga_R_gv.png')
        plt.show() 



col = ['orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue',\
'red','green','orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue',\
'red','green']

if __name__ == "__main__":
    main()
