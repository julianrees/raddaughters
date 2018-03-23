library(reshape2)
library(plyr)
library(drc)
library(drfit)
library(ggplot2)
library(deSolve)



#Specific activity constants
activity = c(k1 = 58, k2 = 180000, k3 = 1600000000, k4 = 20000, k5 = 1.3*10^13, k6 = 4700, k7 = 0, k8 = 96216216216, k9 = 410000) #these are specific activities in uCi/ng

#initial activity of actinium-225 in uCi
uciac225 = 0.0001 #thats 0.1 nCi = 220 CPM

#Initial nmoles of actinium-225
nmolesac225 = uciac225/58/225

#initial DPM of ac-225
dpmac225 = uciac225*2200000

#lambda values are ln(2)/t(1/2) with time in days since *24, other wise in hours
parameters = c(l1 = 0.002888113*24, l2 = 8.487516497*24, l3 = 77254.79412*24, l4 = 0.904105018*24, l5 = 594126154.8*24, l6 = 0.213276056*24, l7 = 3.64*10^-20*24, l8 = 4620981.204*24, l9 = 19.24517854*24)

#probabilities of the destruction to a certain species
probabilities = c(ac2fr = 1, fr2at = 1, at2bi = 0.99923, bi2po = 0.978, po2pb = 1, pb2bi = 1, at2rn = 0.00077, bi2tl = 0.0209, rn2po = 1, tl2pb = 1)

#these numbers are in nmoles
state = c(A = nmolesac225, B = 0, C = 0, D = 0, E = 0, f = 0, G = 0, H = 0, I = 0)

masses = c(j1 = 225, j2 = 221, j3 = 217, j4 = 213, j5 = 213, j6 = 209, j7 = 209, j8 = 217, j9 = 209)


#calculate nmoles of species
daughters = function(t, state, parameters, probabilities) {with(as.list(c(state, parameters, probabilities)),{
  dA = -l1*A
  dB = l1*ac2fr*A-l2*B
  dC = l2*fr2at*B-l3*C
  dH = l3*at2rn*C-l8*H
  dD = l3*at2bi*C-l4*D
  dI = l4*bi2tl*D-l9*I
  dE = l4*bi2po*D+l8*rn2po*H-l5*E
  df = l5*po2pb*E+l9*tl2pb*I-l6*f
  dG = l6*pb2bi*f-l7*G
  
  
  list(c(dA, dB, dC, dD, dE, df, dG, dH, dI))
})}

timedays = 5
timestep = 0.001
timestepout = 1/timestep

times = seq(0, timedays, by = timestep)
timesout = seq(1, timedays*timestepout+1, by = 1)

#MODEL INTEGRATION

out = ode(y = state, times = times, func = daughters, parms = parameters, prob=probabilities)
head(out)


#calculate activity produces over time
#daughtersactiv = function(t, activity, masses, out) {with(as.list(c(activity, masses, out)),{

#multiply to get DPM from uCi and divide by initial ac-225 DPM to get Fraction Activity Remaining
Ac225 = masses[1]*activity[1]*out[timesout,2]*2200000/dpmac225
Fr221 = masses[2]*activity[2]*out[timesout,3]*2200000/dpmac225
At217 = masses[3]*activity[3]*out[timesout,4]*2200000/dpmac225
Bi213 = masses[4]*activity[4]*out[timesout,5]*2200000/dpmac225
Po213 = masses[5]*activity[5]*out[timesout,6]*2200000/dpmac225
Pb209 = masses[6]*activity[6]*out[timesout,7]*2200000/dpmac225
Bi209 = masses[7]*activity[7]*out[timesout,8]*2200000/dpmac225
Rn217 = masses[8]*activity[8]*out[timesout,9]*2200000/dpmac225
Tl209 = masses[9]*activity[9]*out[timesout,10]*2200000/dpmac225
SUM = (Ac225+Fr221+At217+Bi213+Po213+Pb209+Bi209+Rn217+Tl209)
SUMoverac225 = SUM/Ac225

daughtersdata = data.frame(times)
daughtersdata = cbind(daughtersdata, Ac225, Fr221, At217, Bi213, Po213, Pb209, Bi209, Rn217, Tl209, SUM, SUMoverac225)
colnames(daughtersdata)[12] = "SUM / Ac225"


#melt this first
mdaughtersdata = melt(daughtersdata, id="times")

#plot the indivudual activities produced

ggplot(mdaughtersdata, aes(x=times, y=value, by=variable))+
    geom_point(aes(color=variable), size = 1.25)+
  
  
  scale_x_log10(breaks=c(0.0001, 0.001, 0.01, 0.1, 1, 10))+
  annotation_logticks(base = 10, sides = "b", scaled = TRUE,
                      short = unit(0.1, "cm"), mid = unit(0.2, "cm"), long = unit(0.3, "cm"),
                      colour = "black", size = 0.5, linetype = 1, alpha = 1, color = NULL)+
  
  scale_y_continuous(labels = scales::percent, breaks=c(1, 2, 3, 4, 5, 6))+

  
  labs(x = "Time (days)", y = "% Activity(t) / Ac225(0)")+
  theme(text = element_text(size=18))+
  guides(color=guide_legend(title=""))
  


#melt the masses
out = data.frame(out)
colnames(out) = c("time", "Ac-225", "Fr-221", "At-217", "Bi-213", "Po-213", "Pb-209", "Bi-209", "Rn-217", "Tl-209")
mout = melt(out, id="time") 

#plot the masses

ggplot(mout, aes(x=time, y=value, by=variable))+
  geom_point(aes(color=variable))+
  scale_x_log10(breaks=c(0.001, 0.01, 0.1, 1, 10))+
  annotation_logticks(base = 10, sides = "b", scaled = TRUE,
                      short = unit(0.1, "cm"), mid = unit(0.2, "cm"), long = unit(0.3, "cm"),
                      colour = "black", size = 0.5, linetype = 1, alpha = 1, color = NULL)+
  scale_y_log10(breaks=c(10^(-20), 10^(-15), 10^(-10), 10^(-5)))+
  #annotation_logticks(base = 10, sides = "l", scaled = TRUE,
                      #short = unit(0.1, "cm"), mid = unit(0.2, "cm"), long = unit(0.3, "cm"),
                      #colour = "black", size = 0.5, linetype = 1, alpha = 1, color = NULL)+
  labs(x = "Time (days)", y = "Amount (nmoles)")+
  theme(text = element_text(size=15))+
  guides(color=guide_legend(title=""))



#what if start out with actinium, do TLC plate, and now the daughters are separated? Start initial amounts @ time = 2/24, or 2 hours after blotting

#
##
####
#######
###############
################################
###################################################################################################################################################
################################
###############
#######
####
##
#


#These values for 'state' are using 'out' from the above calcualtion, so need to run everything!
parameters = c(l1 = 0.002888113*24, l2 = 8.487516497*24, l3 = 77254.79412*24, l4 = 0.904105018*24, l5 = 594126154.8*24, l6 = 0.213276056*24, l7 = 3.64*10^-20*24, l8 = 4620981.204*24, l9 = 19.24517854*24)
state = c(A = 0, B = out[out$time==0.083,"Fr-221"], C = out[out$time==0.083,"At-217"], D = out[out$time==0.083,"Bi-213"], E = out[out$time==0.083,"Po-213"], f = out[out$time==0.083,"Pb-209"], G = out[out$time==0.083,"Bi-209"], H = out[out$time==0.083,"Rn-217"], I = out[out$time==0.083,"Tl-209"]) #these numbers are in nmoles
activity = c(k1 = 58, k2 = 180000, k3 = 1600000000, k4 = 20000, k5 = 1.3*10^13, k6 = 4700, k7 = 0, k8 = 96216216216, k9 = 410000)
masses = c(j1 = 225, j2 = 221, j3 = 217, j4 = 213, j5 = 213, j6 = 209, j7 = 209, j8 = 217, j9 = 209)
#initmassng = c(i1 = 1, i2 = 0,  i3 = 0,  i4 = 0,  i5 = 0,  i6 = 0,  i7 = 0)

#calculate nmoles of species
daughters = function(t, state, parameters, probabilities) {with(as.list(c(state, parameters, probabilities)),{
  dA = -l1*A
  dB = l1*ac2fr*A-l2*B
  dC = l2*fr2at*B-l3*C
  dH = l3*at2rn*C-l8*H
  dD = l3*at2bi*C-l4*D
  dI = l4*bi2tl*D-l9*I
  dE = l4*bi2po*D+l8*rn2po*H-l5*E
  df = l5*po2pb*E+l9*tl2pb*I-l6*f
  dG = l6*pb2bi*f-l7*G
  
  
  list(c(dA, dB, dC, dD, dE, df, dG, dH, dI))
})}

#timedays = 10
#timestep = 0.001
#timestepout = 1/timestep

#times = seq(0, timedays, by = timestep)
#timesout = seq(1, timedays*timestepout+1, by = 1)

#MODEL INTEGRATION

out = ode(y = state, times = times, func = daughters, parms = parameters, prob=probabilities)
head(out)


#calculate activity produces over time
#daughtersactiv = function(t, activity, masses, out) {with(as.list(c(activity, masses, out)),{

Ac225 = masses[1]*activity[1]*out[timesout,2]*2200000/dpmac225
Fr221 = masses[2]*activity[2]*out[timesout,3]*2200000/dpmac225
At217 = masses[3]*activity[3]*out[timesout,4]*2200000/dpmac225
Bi213 = masses[4]*activity[4]*out[timesout,5]*2200000/dpmac225
Po213 = masses[5]*activity[5]*out[timesout,6]*2200000/dpmac225
Pb209 = masses[6]*activity[6]*out[timesout,7]*2200000/dpmac225
Bi209 = masses[7]*activity[7]*out[timesout,8]*2200000/dpmac225
Rn217 = masses[8]*activity[8]*out[timesout,9]*2200000/dpmac225
Tl209 = masses[9]*activity[9]*out[timesout,10]*2200000/dpmac225
SUM = (Ac225+Fr221+At217+Bi213+Po213+Pb209+Bi209+Rn217+Tl209)
#SUMoverac225 = SUM/Ac225

daughtersdata = data.frame(times)
daughtersdata = cbind(daughtersdata, Ac225, Fr221, At217, Bi213, Po213, Pb209, Bi209, Rn217, Tl209, SUM)#, SUMoverac225)
colnames(daughtersdata) = c("times", "Ac-225", "Fr-221", "At-217", "Bi-213", "Po-213", "Pb-209", "Bi-209", "Rn-217", "Tl-209", "SUM")



daughtersdata <- daughtersdata[c(1,2,3,4,5,6,7,9,10,8,11)]


#colnames(daughtersdata)[12] = "SUM / Ac225"

#melt this first
mdaughtersdata = melt(daughtersdata, id="times")
colnames(mdaughtersdata) <- c("times","Species","value")



#plot the indivudual activities produced



ggplot(mdaughtersdata, aes(x=times, y=value, by=Species))+
  geom_point(aes(color=Species, shape=Species), size=1.25, alpha=1, stroke = 1.25)+
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17))+ 
  scale_x_log10(breaks=c(0.001, 0.01, 0.1, 1, 10))+
  annotation_logticks(base = 10, sides = "bl", scaled = TRUE,
                      short = unit(0.1, "cm"), mid = unit(0.2, "cm"), long = unit(0.3, "cm"),
                      colour = "black", size = 0.5, linetype = 1, alpha = 1, color = NULL)+
  
  scale_y_log10(labels = scales::percent, breaks=c(10^(-4):1 %o% 10^(1:4)), limits = c(2*10^(-4),1))+
  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  
  labs(x = "Time (days)", y = "% Activity(t) / Ac-225(0.1), w/o Ac-225 Present",  color="Species")+
  theme(text = element_text(size=18, face = "bold"),
        axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"))+
  guides(shape=guide_legend(override.aes = list(size=3)))
  



  


#melt the masses
out = data.frame(out)
colnames(out) = c("time", "Ac225", "Fr221", "At217", "Bi213", "Po213", "Pb209", "Bi209", "Rn217", "Tl209")
mout = melt(out, id="time") 

#plot the masses

ggplot(mout, aes(x=time, y=value, by=variable))+
  geom_point(aes(color=variable))+
  scale_x_log10(breaks=c(0.001, 0.01, 0.1, 1, 10))+
  annotation_logticks(base = 10, sides = "b", scaled = TRUE,
                      short = unit(0.1, "cm"), mid = unit(0.2, "cm"), long = unit(0.3, "cm"),
                      colour = "black", size = 0.5, linetype = 1, alpha = 1, color = NULL)+
  scale_y_log10(breaks=c(10^(-15), 10^(-10), 10^(-5), 10^(-1)))+
  #annotation_logticks(base = 10, sides = "l", scaled = TRUE,
  #short = unit(0.1, "cm"), mid = unit(0.2, "cm"), long = unit(0.3, "cm"),
  #colour = "black", size = 0.5, linetype = 1, alpha = 1, color = NULL)+
  labs(x = "Time (days)", y = "Amount (nmoles)")+
  theme(text = element_text(size=15))+
  guides(color=guide_legend(title=""))


#### Activity determination of 229Th by means of liquid scintillation counting - 2013 - Kossert ####
