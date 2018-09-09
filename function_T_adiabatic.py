# -*- coding: utf-8 -*-
"""
Determine adibatic temperature combustion (with and without dissociation)
@author: Kevin Sartor
"""
#Libraries import
from Ideal_gaz_properties import IG_props,molarMass
from scipy.optimize import fsolve
import math
#import scipy
import numpy as np

def T_adiab(e):
        
    n=16
    X = np.zeros(n)
    
    "Mass dry composition"
    Cm=50.2                             #Mass fraction of C"
    Hm =6.1                             #Mass fraction of H"
    Nm=0.09*1                           #Mass fraction of N"
    Sm=0.025                           #Mass fraction of S"
    Ashm=0.3                          #Mass fraction of Ash"
    Om = 100-(Cm+Hm+Nm+Sm+Ashm)         #Mass fraction of O (usually calculated by difference EN15104"
    PCI_boie_mas=17432849               #PCI massique du combustible
    PCI_boie_mas=17624469.502914
    
    "Mass wet composition"
    H2O=7.6*1                             #Mass fraction of H2O"
    C = Cm*(100-H2O)/100                # 
    H = Hm*(100-H2O)/100                # 
    O = Om*(100-H2O)/100                # 
    N = Nm*(100-H2O)/100                # 
    S = Sm*(100-H2O)/100                #
    Ash=Ashm*(100-H2O)/100              # 
    #PCS_boie_mas_emp=1E6*((34.91*C/100+117.83*H/100-10.34*O/100+10.05*S/100-1.51*N/100)*(1-Ash/100)-2.11*Ash/100)/(1-Ash/100)
    #PCI_boie_mas_emp=35160*Cm/100+94438*Hm/100-11090*Om/100+6280*Nm/100+10465*Sm/100
    #print (PCS_boie_mas_emp)
    "Mole fraction"
    m = C/(100*molarMass['C'])
    n = H/(100*molarMass['H'])
    x = O/(100*molarMass['O'])
    y = N/(100*molarMass['N'])
    z = S/(100*molarMass['S'])
    x_water=H2O/(100*molarMass['H2O'])
    
    "Molar fraction of component in the fuel"
    Y_m=m/(m+n+x+y+z+x_water)
    Y_n=n/(m+n+x+y+z+x_water)
    Y_x=x/(m+n+x+y+z+x_water)
    Y_y=y/(m+n+x+y+z+x_water)
    Y_z=z/(m+n+x+y+z+x_water)
    Y_water=x_water/(m+n+x+y+z+x_water)
    #print (Y_m,Y_n,Y_x,Y_y,Y_z,Y_water)
    "Air composition"
    T_air=25+273.15	        
    rh=.5
    #P_comb =101325	#[Pa] "Combustion pressure = the same as the atmospheric pressure"
    #P_atm=101325
    WAR=.009882#HumRat(AirH2O,T=T_air,r=rh,P=P_comb)
    pc_O2=.20948
    pc_Ar=0.009*0
    pc_N2=1-pc_Ar-pc_O2
    
    # Mass fraction of dry air
    x_bar_air_dry_O2=pc_O2*molarMass['O2'] /(pc_O2*molarMass['O2']+pc_N2*molarMass['N2']+pc_Ar*molarMass['Ar'])
    x_bar_air_dry_N2=pc_N2*molarMass['N2'] /(pc_O2*molarMass['O2']+pc_N2*molarMass['N2']+pc_Ar*molarMass['Ar'])
    x_bar_air_dry_Ar=pc_Ar*molarMass['Ar'] /(pc_O2*molarMass['O2']+pc_N2*molarMass['N2']+pc_Ar*molarMass['Ar'])
    
    # Mass fraction of wet air
    x_bar_air_wet_O2=x_bar_air_dry_O2/(1+WAR);
    x_bar_air_wet_N2=x_bar_air_dry_N2/(1+WAR);
    x_bar_air_wet_Ar=x_bar_air_dry_Ar/(1+WAR);
    x_bar_air_wet_H2O=WAR/(1+WAR)
    
    # Molar fraction of wet air
    x_air_wet_H2O=x_bar_air_wet_H2O/molarMass['H2O']/(x_bar_air_wet_O2/molarMass['O2']+x_bar_air_wet_N2/molarMass['N2']+x_bar_air_wet_H2O/molarMass['H2O']+x_bar_air_wet_Ar/molarMass['Ar'])
    x_air_wet_N2=x_bar_air_wet_N2/molarMass['N2']/(x_bar_air_wet_O2/molarMass['O2']+x_bar_air_wet_N2/molarMass['N2']+x_bar_air_wet_H2O/molarMass['H2O']+x_bar_air_wet_Ar/molarMass['Ar'])
    x_air_wet_O2=x_bar_air_wet_O2/molarMass['O2']/(x_bar_air_wet_O2/molarMass['O2']+x_bar_air_wet_N2/molarMass['N2']+x_bar_air_wet_H2O/molarMass['H2O']+x_bar_air_wet_Ar/molarMass['Ar'])
    x_air_wet_Ar=x_bar_air_wet_Ar/molarMass['O2']/(x_bar_air_wet_O2/molarMass['O2']+x_bar_air_wet_N2/molarMass['N2']+x_bar_air_wet_H2O/molarMass['H2O']+x_bar_air_wet_Ar/molarMass['Ar'])
    
    "Combustion parameters"
#    e=.015                               # Air excess
    efficiency=1                        # Combustion efficiency (should be calculated by CO, CH4, NO losses)
    P=101325                            # Combustion pressure
    lambda_a=1+e
    #phi=1/lambda_a
    #f_st =1/((m+n/4-x/2+z)*(molarMass['O2']+(x_air_wet_N2/x_air_wet_O2)*molarMass['N2']+x_air_wet_H2O/x_air_wet_O2*molarMass['H2O']))
    #f = f_st/(1+e)
    #AFR=1/f
    
    "Molar Enthalpy reference"
    h_CO2_ref = IG_props('H','CO2',298.15,1)
    h_CO_ref  = IG_props('H','CO',298.15,1)
    h_O2_ref  = IG_props('H','O2',298.15,1)
    h_H2O_ref = IG_props('H','H2O',298.15,1)
    h_H2_ref  = IG_props('H','H2',298.15,1)
    h_NO_ref  = IG_props('H','NO',298.15,1)
    h_N2_ref  = IG_props('H','N2',298.15,1)
    h_OH_ref  = IG_props('H','OH',298.15,1)
    h_H_ref   = IG_props('H','H',298.15,1)
    h_O_ref   = IG_props('H','O',298.15,1)
    h_CH4_ref = IG_props('H','CH4',298.15,1)
    h_N_ref   = IG_props('H','N',298.15,1)
    h_NO2_ref = IG_props('H','NO2',298.15,1)
    h_SO2_ref = IG_props('H','SO2',298.15,1)
    h_SO3_ref = IG_props('H','SO3',298.15,1)
    h_Ar_ref  = IG_props('H','Ar',298.15,1)
    
    "Molar Enthalpy of reactives"
    h_CO2_r = IG_props('H','CO2',T_air,1)
    h_CO_r  = IG_props('H','CO',T_air,1)
    h_O2_r  = IG_props('H','O2',T_air,1)
    h_H2O_r = IG_props('H','H2O',T_air,1)
    h_H2_r  = IG_props('H','H2',T_air,1)
    h_NO_r  = IG_props('H','NO',T_air,1)
    h_N2_r  = IG_props('H','N2',T_air,1)
    h_OH_r  = IG_props('H','OH',T_air,1)
    h_H_r   = IG_props('H','H',T_air,1)
    h_O_r   = IG_props('H','O',T_air,1)
    h_CH4_r = IG_props('H','CH4',T_air,1)
    h_N_r   = IG_props('H','N',T_air,1)
    h_NO2_r = IG_props('H','NO2',T_air,1)
    h_SO2_r = IG_props('H','SO2',T_air,1)
    h_SO3_r = IG_props('H','SO3',T_air,1)
    h_Ar_r  = IG_props('H','Ar',T_air,1)
    
    "Calcul PCI Molaire"
    M_fuel=Y_m*molarMass['C']+Y_n*molarMass['H']+Y_x*molarMass['O']+Y_y*molarMass['N']+Y_z*molarMass['S']+Y_water*molarMass['H2O']
    PCI_boie_mol=PCI_boie_mas*M_fuel
    #PCI_boie_mol=139918109.55
    
    "Calcul des coefficients "
    a_s=1                                                              # 1 mole de combustible
    b_s=a_s*Y_m                                                        # equilibre C
    j_s=a_s*Y_z                                                        # equilibre sur S
    n_s=Y_m+Y_n/4-Y_x/2+Y_z                                            # coef O2
    d_s=n_s*e                                                          # equilibre sur O2
    e_s=n_s*(1+e)*(x_air_wet_N2/x_air_wet_O2)+a_s*Y_y/2                # equilibre sur N2
    c_s=(a_s*Y_n/2+Y_water+n_s*(1+e)*(x_air_wet_H2O/x_air_wet_O2)) # equilibre sur H2O
    #print (b_s,c_s,d_s,e_s,n_s,j_s)
    
    "Adiabatic combustion"
    def combustion_wo_diss(p): 
#        global a_s,b_s,c_s,d_s,e_s,n_s,j_s
        T  = p 
        h_CO2_T = IG_props('H','CO2',T,1)
        h_O2_T  = IG_props('H','O2',T,1)
        h_H2O_T = IG_props('H','H2O',T,1)
        h_N2_T  = IG_props('H','N2',T,1)
        h_Ar_T  = IG_props('H','Ar',T,1)
        H_r=(a_s*(PCI_boie_mol*efficiency)+n_s*(1+e)*(h_O2_r-h_O2_ref)+n_s*(1+e)*((x_air_wet_N2/x_air_wet_O2))*(h_N2_r-h_N2_ref)+((x_air_wet_H2O/x_air_wet_O2))*(h_H2O_r-h_H2O_ref)+n_s*(1+e)*(x_air_wet_Ar/x_air_wet_O2)*(h_Ar_r-h_Ar_ref))
        residu = b_s*(h_CO2_T-h_CO2_ref)+c_s*(h_H2O_T-h_H2O_ref)+e_s*(h_N2_T-h_N2_ref)+d_s*(h_O2_T-h_O2_ref)+n_s*(1+e)*(x_air_wet_Ar/x_air_wet_O2)*(h_Ar_T-h_Ar_ref)-H_r#(a_s*(PCI_boie_mol*efficiency)+n_s*(1+e)*(h_O2_r-h_O2_ref)+n_s*(1+e)*((x_air_wet_N2/x_air_wet_O2))*(h_N2_r-h_N2_ref)+((x_air_wet_H2O/x_air_wet_O2))*(h_H2O_r-h_H2O_ref)+n_s*(1+e)*(x_air_wet_Ar/x_air_wet_O2)*(h_Ar_r-h_Ar_ref))
        return (residu)
    Temp_adiab = fsolve(combustion_wo_diss,(1500))                      # Guess value = 1500
    n_r=a_s+n_s*(1+e)+n_s*(1+e)*(x_air_wet_N2/x_air_wet_O2)+(x_air_wet_H2O/x_air_wet_O2)*n_s*(1+e)+n_s*(1+e)*(x_air_wet_Ar/x_air_wet_O2)
    n_p=b_s+c_s+d_s+e_s+j_s+n_s*(1+e)*(x_air_wet_Ar/x_air_wet_O2)
    contrib_wet_air=(1+e)*n_s*x_air_wet_H2O/x_air_wet_O2
    #print (contrib_wet_air)
    n_tot=(Y_n+2*Y_water+2*contrib_wet_air)/(2*c_s/n_p)
    
    #n=nn/mm
    #Y_m=mm/(nn+mm)
    #Y_n=mm/(nn+mm)
    #mm=pc_C/(100*12)
    #nn=pc_H/(100*1)
    #pc_H+pc_C=100
    
    T_s=Temp_adiab[0]
    #n_t=b_s+c_s+d_s+e_s+n_s*(1+e)*(x_air_wet_Ar/x_air_wet_O2)
    #print (n_p,n_r)
#    print ("T sans dissociation: "+str(round(T_s,1))+" K")
    T=T_s
    T_diss=T_s-15
    R=8314
    #boucle de convergence on cherche à minimiser le DT mais on limite à 50 itérations.
    boucle=0
    while (math.fabs(T_diss-T)>0.1 and boucle<50):
        T=T_diss
        contrib_wet_air=(1+e)*(m+n/4-x/2+z)*x_air_wet_H2O/x_air_wet_O2
        RAO = (x_air_wet_Ar)/(contrib_wet_air+x+x_water+(m+n/4-x/2+z)*2*(1+e))
        RCO = (m)/(contrib_wet_air+x_water+x+(m+n/4-x/2+z)*2*(1+e))
        RHO = (2*contrib_wet_air+n+2*x_water)/(contrib_wet_air+x+x_water+(m+n/4-x/2+z)*2*(1+e))
        RNO = (y+lambda_a*(m+n/4-x/2+z)*2*(x_air_wet_N2)/x_air_wet_O2)/(contrib_wet_air+x+x_water+(m+n/4-x/2+z)*2*(1+e))
        
        PT=P/101325.
        RNC=RNO/RCO
        ROC=1./RCO
        RHC=RHO/RCO
        RAC=RAO/RCO
        
        CONO=1./(RCO+RHO+RNO+RAO+1)
        CONN=CONO*RNO
        CONA=CONO*RAO
        CONC=CONO*RCO
        CONH=CONO*RHO
        PHI=2*RCO+RHO/2
        if PHI==1:
            PHI=1.00001
        RHTH2O=math.tanh(10**(0.002353*(5000.0-T)-4.0))
        RHTCO2=math.tanh(10**(-0.00126262626*(T-4417.0)-2.0))
        if T<=2600:
            RHTCO2=1.0
            RHTH2O=1.0
        DIV=4.4375*PHI**2-8.875*PHI+7.0975    
        
        PEAK=0.1347*PHI+0.00459
        QUAL=(1-PHI+1e-9)/math.fabs(1-PHI+1e-9)
        CO2INT=(1+QUAL)*PEAK/2+(1.0-QUAL)*(0.28-PEAK)/2
        MULT=1
        if(PHI>1.8):
            MULT=(2.0-PHI)/0.2
        PPCO2=((math.sqrt(PT*CO2INT*math.tanh(1/RHC))/DIV)*MULT+(PT*0.2/((1/RHC)+RHC))*(1-MULT))*RHTCO2
        PPH2O=PT*0.2*(math.tanh(RHC/10))*RHTH2O
        
        CNT=CONO/2+CONH/4+CONN/2+CONA
        PPO2=PT*(CONO/2-CONC-CONC/4)/CNT
        
        g_CO2   = IG_props('G','CO2',T,1)
        g_CO    = IG_props('G','CO',T,1)
        g_O2    = IG_props('G','O2',T,1)
        g_H2O   = IG_props('G','H2O',T,1)
        g_H2    = IG_props('G','H2',T,1)
        g_NO    = IG_props('G','NO',T,1)
        g_N2    = IG_props('G','N2',T,1)
        g_OH    = IG_props('G','OH',T,1)
        g_H     = IG_props('G','H',T,1)
        g_O     = IG_props('G','O',T,1)
        g_CH4   = IG_props('G','CH4',T,1)
        g_N     = IG_props('G','N',T,1)
        g_NO2   = IG_props('G','NO2',T,1)
        g_SO2   = IG_props('G','SO2',T,1)
        g_SO3   = IG_props('G','SO3',T,1)
        
        BK1=math.exp(-((1/2*g_O2+g_CO)-g_CO2)/(R*T))
        BK2=math.exp(-((1/2*g_O2+g_H2)-g_H2O)/(R*T))
        BK3=math.exp(-((g_OH+1/2*g_H2)-(g_H2O))/(R*T))
        BK4=math.exp(-((g_H)-(1/2*g_H2))/(R*T))
        BK5=math.exp(-((g_O)-(1/2*g_O2))/(R*T))
        BK7=math.exp(-((g_CH4+2*g_H2O)-(g_CO2+4*g_H2))/(R*T))
        BK9=math.exp(-((g_N)-(1/2*g_N2))/(R*T))
        BK6=math.exp(-(g_NO-1/2*g_N2-1/2*g_O2)/(R*T))
        BK8=math.exp(-((g_NO2)-(1/2*g_N2+g_O2))/(R*T))
        
        #bk=10**(-2.602+2*5.725-13.801)
        #print (bk,BK7)
        #Start iteration
        DELH2O=0
        DELH2=0
        DELO2=0
        DELCO2=0
        I=1
        stop=0
        while (I<50 and stop==0):
            JA=I
            KA=0
           
            PP1CO2=PPCO2
            PP1H2O=PPH2O
            PPH2O=PPH2O+DELH2O
            if (PPH2O<0):
                PPH2O=PP1H2O/2
                KA=1
            #4
            PPCO2=PPCO2+DELCO2
            if(PPCO2<0):
               PPCO2=PP1CO2/2
               KA=1 
           
            if (PHI>=1.0):
                print ('Riche')
                PP1H2=PPH2
                PPH2=PPH2+DELH2
                if(PPH2<0):
                    PPH2=PP1H2/2.0
                    KA=1
                    PPO2=(PPH2O*BK2/PPH2)**2.0
                    PPCH4=BK7*PPCO2*(PPH2**4)/(PPH2O*PPH2O)
            else:
                PP1O2=PPO2
                PPO2=PPO2+DELO2
                if I<2:
    #                print ('A')#(PPO2)
                    tempo=1
                if (PPO2<=0):
                    PPO2=PP1O2/2.
                    KA=1
                PPH2=PPH2O*BK2/math.sqrt(PPO2)
                PPCH4=0.0
        
            PPCO=PPCO2*BK1/math.sqrt(PPO2)
            PPOH=PPH2O*BK3/math.sqrt(PPH2)
            PPO=BK5*math.sqrt(PPO2)
            PPH=BK4*math.sqrt(PPH2)
            PPA=RAC*(PPCO2+PPCO+PPCH4)
            PPNO=BK6**2*PPO2*(math.sqrt(1.0+8.0*RNC*(PPCO2+PPCO+PPCH4)/(PPO2*BK6**2))-1.0)/4.0
            PPN2=0.5*RNC*(PPCO2+PPCO+PPCH4)-0.5*PPNO
            PPN=BK9*PPN2**0.5
            PPNO2=BK8*PPO2*PPN2**(0.5)
            PA=PT-PPNO2-PPN 
            
            C1=2.0+(1.0/PPH2O)*(2.0*PPH2+.5*PPOH+.5*PPH)
            C2=(1.0/PPO2)*(-PPH2+.25*PPOH-.25*PPH+.5*RHC*PPCO)
            C3=-RHC*(1.0+PPCO/PPCO2)
            C4=2.0*PPH2O+2.0*PPH2+PPOH+PPH-RHC*(PPCO2+PPCO)
            C5=1.0+.5*PPOH/PPH2O
            C6=2.0+(-.5*PPCO+.25*PPOH+.5*PPO+.5*ROC*PPCO)/ PPO2+(PPNO/PPO2) * (2.0*PPN2-.5*RNC*PPCO)/(4.0*PPN2+PPNO)
            C7=2.0-ROC-ROC*(PPCO/PPCO2)+(PPCO/PPCO2)+(RNC*PPNO/PPCO2)*(PPCO+PPCO2)/(4.0*PPN2+PPNO)
            C8=2.0*PPCO2+PPCO+PPH2O+2.0*PPO2+PPOH+PPO-ROC*(PPCO2+PPCO)+PPNO
            C9=1.0+(1.0/PPH2O)*(PPH2+.5*PPOH+.5*PPH)+RAO*(1.0+PPOH/PPH2O)
            C10=1.0+(-0.5*PPCO-.5*PPH2+.25*PPOH-.25*PPH+.5*PPO-(PPN2*(RNC*PPCO+PPNO)+PPNO*(2.0*PPN2-.5*RNC*PPCO))/(4.0*PPN2+PPNO))/PPO2-0.5*RAC*PPCO/PPO2
            C11=1.0+(PPCO+2.0*RNC*(PPN2*(PPCO2+PPCO)+.5*PPNO*(PPCO2+PPCO)) /(4.0*PPN2+PPNO))/PPCO2+RAC*(1.0+PPCO/PPCO2)
            C12=PPCO2+PPCO+PPH2O+PPO2+PPH2+PPOH+PPH+PPO+PPN2+PPNO+PPA-PA
            
            DET1=C2*(C7*C12-C11*C8)-C3*(C6*C12-C10*C8)+C4*(C6*C11-C1*C10)
            DET2=C1*(C7*C12-C11*C8)-C3*(C5*C12-C9 *C8)+C4*(C5*C11-C9*C7 )
            DET3=C1*(C6*C12-C10*C8)-C2*(C5*C12-C9 *C8)+C4*(C5*C10-C6*C9 )
            DET4=C1*(C6*C11-C10*C7)-C2*(C5*C11-C9 *C7)+C3*(C5*C10-C9*C6 )
            
            if (DET4==0):
                DET4=1.0E-100
            DELCO2=-DET3/DET4
            
            if(PHI<1.0):
                DELH2O=-DET1/DET4
                DELO2=DET2/DET4
            #    IF(K==1}GOTO 14
            #    IF(I<3)GOTO 14
            #    IF(ABS(DEL02/PP02)-1.E-6)13,13,14
            #12 CONTINU
            else:
                DELH2=-DET1/DET4
                DELH2O=DET2/DET4
                print ('Error')
            #if(KA==1)GOTO 14
        #    if(I>=3):
                #IF(ABS(DELH2/PPH2).GT.1.E-6)GOTO 14
            #13
            #IF(ABS(DELH20/PPB20).GT.1.E-6)GOTO 14
            #IF(ABS(DELC02/ PPC02).GT.1.E-6)GOTO 14
        #    time.sleep(0.5)
            if(math.fabs(C12/PA) < 1.E-6):
                stop=1
            I=I+1
            #14 CONTINUE
            #PRINT 120
            #15
        
            X[1]=PPH2/PT
            X[2]=PPO2/PT
            X[3]=PPH2O/PT
            X[4]=PPCO/PT
            X[5]=PPCO2/PT
            X[6]=PPOH/PT
            X[7]=PPH/PT
            X[8]=PPO/PT
            X[9]=PPN2/PT
            X[10]=PPN/PT
            X[11]=PPNO/PT
            X[12]=PPNO2/PT
            X[13]=PPCH4/PT
            X[14]=PPA/PT
        
    #    print (np.sum(X))
        n_tot=(Y_n+2*Y_water+2*(1+e)*n_s*x_air_wet_H2O/x_air_wet_O2)/((2* X[1]+2* X[3]+ X[6]+ X[7]+4* X[13]))
    #    print (n_tot)
        n_p=n_tot
        
        b_d=X[5]*n_p#*(p_tot/p_0)
        f_d=X[4]*n_p#-(f_d/n_p)#*(p_tot/p_0)
        d_d=X[2]*n_p#-(d_d/n_p)#*(p_tot/p_0)
        c_d=X[3]*n_p#-(c_d/n_p)#*(p_tot/p_0)
        i_d=X[1]*n_p#-(i_d/n_p)#*(p_tot/p_0)
        h_d=X[6]*n_p#-(h_d/n_p)#*(p_tot/p_0)
        k_d=X[7]*n_p#-(k_d/n_p)#*(p_tot/p_0)
        l_d=X[8]*n_p#-(l_d/n_p)#*(p_tot/p_0)
        e_d=X[9]*n_p#-(e_d/n_p)#*(p_tot/p_0)
        q_d=X[10]*n_p#-(q_d/n_p)#*(p_tot/p_0)
        g_d=X[11]*n_p#-(g_d/n_p)#*(p_tot/p_0)
        o_d=X[12]*n_p#-(o_d/n_p)#*(p_tot/p_0)
        m_d=X[13]*n_p#-(m_d/n_p)#*(p_tot/p_0)
        j_d=0
        p_d=0
        
        LHV_CO=-(h_CO2_ref-h_CO_ref-1/2*h_O2_ref)
        LHV_NO=-(h_NO_ref-1/2*h_N2_ref-1/2*h_O2_ref)
        LHV_CH4=-((h_CH4_ref+2*h_H2O_ref)-(h_CO2_ref+4*h_H2_ref))
        LHV_H2=h_H2_ref+1/2*h_O2_ref-h_H2O_ref
        LHV_NO2=-((h_NO2_ref)-(1/2*h_N2_ref+h_O2_ref))
        q_loss_CO=f_d*LHV_CO
        q_loss_NO=g_d*LHV_NO
        q_loss_NO2=o_d*LHV_NO2
        q_loss_CH4=m_d*LHV_CH4
        q_loss_H2=i_d*LHV_H2
        q_loss=q_loss_CO+q_loss_CH4+q_loss_H2
        q_loss_CO_mas=q_loss_CO/28
        q_loss_H2_mas=q_loss_H2/2
        q_loss_mas=q_loss_CO_mas+q_loss_H2_mas
        q_loss_mas=q_loss_CO_mas+q_loss_H2_mas
        epsilon_comb=(PCI_boie_mas-n_p*q_loss_mas)/(PCI_boie_mas)
        "Enthalpies"
        def combustion_diss(p): 
#            global a_d,b_d,c_d,d_d,e_d,f_d,g_d,h_d,i_d,j_d,k_d,l_d,m_d,o_d,p_d,q_d
            T=p
            h_CO2_T = IG_props('H','CO2',T,1)
            h_CO_T  = IG_props('H','CO',T,1)
            h_O2_T  = IG_props('H','O2',T,1)
            h_H2O_T = IG_props('H','H2O',T,1)
            h_H2_T  = IG_props('H','H2',T,1)
            h_NO_T  = IG_props('H','NO',T,1)
            h_N2_T  = IG_props('H','N2',T,1)
            h_OH_T  = IG_props('H','OH',T,1)
            h_H_T   = IG_props('H','H',T,1)
            h_O_T   = IG_props('H','O',T,1)
            h_CH4_T = IG_props('H','CH4',T,1)
            h_N_T   = IG_props('H','N',T,1)
            h_NO2_T = IG_props('H','NO2',T,1)
            h_SO2_T = IG_props('H','SO2',T,1)
            h_SO3_T = IG_props('H','SO3',T,1)
            h_Ar_T  = IG_props('H','Ar',T,1)
        
            H_reactifs=a_s*(PCI_boie_mol*efficiency*epsilon_comb)+n_s*(1+e)*(h_O2_r-h_O2_ref)+n_s*(1+e)*((x_air_wet_N2/x_air_wet_O2))*(h_N2_r-h_N2_ref)+((x_air_wet_H2O/x_air_wet_O2))*(h_H2O_r-h_H2O_ref)+n_s*(1+e)*(x_air_wet_Ar/x_air_wet_O2)*(h_Ar_r-h_Ar_ref)
            H_produits=b_d*(h_CO2_T-h_CO2_ref)+c_d*(h_H2O_T-h_H2O_ref)+e_d*(h_N2_T-h_N2_ref)+j_d*(h_SO2_T-h_SO2_ref)+i_d*(h_H2_T-h_H2_ref)+f_d*(h_CO_T-h_CO_ref)+d_d*(h_O2_T-h_O2_ref)+g_d*(h_NO_T-h_NO_ref)+h_d*(h_OH_T-h_OH_ref)+k_d*(h_H_T-h_H_ref)+l_d*(h_O_T-h_O_ref)+m_d*(h_CH4_T-h_CH4_ref)+q_d*(h_N_T-h_N_ref)+o_d*(h_NO2_T-h_NO2_ref)+p_d*(h_SO3_T-h_SO3_ref)+n_s*(1+e)*(x_air_wet_Ar/x_air_wet_O2)*(h_Ar_T-h_Ar_ref) 
        #    print (H_reactifs,H_produits)
            f6=H_reactifs-H_produits
            return (f6) 
        T_diss = fsolve(combustion_diss,(T_s))#,xtol=1e-30)
        boucle=boucle+1
        if boucle>49:
            T_diss=-999
    return T_diss[0],X[5],X[11]
i=1
e=0.1
n=20
A=np.zeros(n)
while i<20:
    print (e,T_adiab(e))
#    A[i]=T_adiab(e)
    i=i+1
    e=e+.1