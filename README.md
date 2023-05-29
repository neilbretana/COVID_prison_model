# COVID_prison_model
# This repository contains the code for the manuscript "Controlling COVID-19 outbreaks in the correctional setting: a mathematical modelling study"

#Running the code
# This code assumes you have C++ and gsl libraries installed in your machine
# To run the code, enter the following script via the console
# gcc -c ./individual.c -lgsl -lgslcblas -lm
# gcc -c -Werror -Wall main.o individual.o -lgsl -lgslcblas -o <instance_name>
# ./<instance_name>.exe ./input_file.xls ./output_file.xls b 

#Input File 
# The input file contains the initial data that is used to populate the model. It also includes the number of days needed for the simulation.

# Opening the input file, the format is as follows:
 <number of simulation days>
 <number of people in the general community>
 <number of inmates in minimum security prisons>  
 <number of prison staff in minimum security prisons>
 <number of healthcare staff in minimum security prisons>
 <number of essential visitors in minimum security prisons>
 <number of family visitors in minimum security prisons>
 <number of inmates in medium security prisons>  
 <number of prison staff in medium security prisons>
 <number of healthcare staff in medium security prisons>
 <number of essential visitors in medium security prisons>
 <number of family visitors in medium security prisons>
 <number of inmates in maximum security prisons>  
 <number of prison staff in maximum security prisons>
 <number of healthcare staff in maximum security prisons>
 <number of essential visitors in maximum security prisons>
 <number of family visitors in maximum security prisons>

# This data can be modified
  
#Baseline model
# The BaselineInmate folder contains the Baseline model. Use the individual.o, main.o, and individual.c files here to run the baseline model as defined below:
# COVID-19 entry via an infected inmate on day 1; inmates can intermingle with other inmates in the same prison area; SARS-CoV-2 delta strain disease-related parameters applied. 
#Variations of this model can be found on the BaselineHS model to set the first COVID-19 infected entrant as a healthcare worker, or BaselinePS to set the first infected entrant as a prison staff.

#Standard PPE model
# The StandardPPE folder contains the Standard PPE model. Use the individual.o, main.o, and individual.c files here to run the baseline model as defined below:
# Standard face masks in use for all inmates, correctional staff, and healthcare staff; This applies a 5% reduction in probability of onward transmission from the source; and a 67% protection of infection for the recipient. This scenario assumes 100% PPE compliance.
  
#PPE + Quarantine + Isolation model
# The PIQincMask folder contains the PPE + Quarantine + Isolation model. Use the individual.o, main.o, and individual.c files here to run the baseline model as defined below:
# Standard face masks are used by inmates everywhere including outside quarantine; new inmates are quarantined for 14 days with PCR tests at day 1 and day 12. If the PCR test returns a positive result (assuming 100% accuracy), the inmate is put into isolation for 14 days; N95 masks are used by staff in the isolation area with an 18% reduction in probability of onward transmission for the source; and an 85% protection from infection for the recipient is applied. This scenario assumes 100% PPE compliance.
  
#Daily RAT model
# The RATdaily folder contains the Daily RAT model. Use the individual.o, main.o, and individual.c files here to run the baseline model as defined below:
# COVID-19 entry assumed to be from 1 infected inmate on day 1; correctional staff and healthcare staff are subjected to RAT testing every day before entering the prison. A pooled RAT sensitivity of 71% and a specificity of 99% was applied. Prison and healthcare staff returning a positive RAT result are assumed to be sent home and subjected to PCR testing within 24 hours.
  
#Standard mask during transit model
# The PPEtransit folder contains the Standard mask during transit model. Use the individual.o, main.o, and individual.c files here to run the baseline model as defined below:
# Standard PPE masks are used by all inmates in-transit with 100% compliance.
  
#Standard mask during transit model
# The PPEtransit folder contains the Standard mask during transit model. Use the individual.o, main.o, and individual.c files here to run the baseline model as defined below:
# Standard PPE masks are used by all inmates in-transit with 100% compliance.
# A variation of this model uses N95 mask during transit and is within the N95 folder of the PPEtransit folder.
  
#Restrict prison transfer model
# The StopMove folder contains the Restrict prison transfer model. Use the individual.o, main.o, and individual.c files here to run the baseline model as defined below:
# Symptomatic inmates are tested using PCR with results returned the next day. If a COVID-19 positive inmate is confirmed, all prison transfers and all court visits are stopped for all inmates across the whole prison system.
  
#Cell isolation model
# The PPELockCell folder contains the Cell isolation model. Use the individual.o, main.o, and individual.c files here to run the baseline model as defined below:
# Symptomatic inmates are isolated in their cell for 14 days with PCR tests at day 1 and day 12; N95 masks are used by staff interacting with isolated inmates; standard PPE masks are used by isolated inmates; isolation occurs on a rolling basis while prison interactions outside isolated cells continue as normal.
  
#Unit isolation model
# The PPELockUnit folder contains the Unit isolation model. Use the individual.o, main.o, and individual.c files here to run the baseline model as defined below:
# Symptomatic inmates are tested using PCR with results by the next day. Units with symptomatic inmates are locked down until no one is actively infected with COVID-19; N95 masks are used by staff interacting with isolated inmates; standard PPE masks are used by isolated inmates; inmates within the unit are free to move within the same unit, but there is no travel to or from another unit within the area of the centre; isolation occurs on a rolling basis while prison interactions outside isolated units continue as normal.
 
#Area isolation model
# The PPELockArea folder contains the Area isolation model. Use the individual.o, main.o, and individual.c files here to run the baseline model as defined below:
# Symptomatic inmates are tested using PCR with results by the next day. Areas with symptomatic inmates are locked down until no one is actively infected with COVID-19; N95 masks are used by staff interacting with isolated inmates; standard PPE masks are used by isolated inmates; inmates within the area are free to move within the area but there is no travel to or from another area within the centre; isolation occurs on a rolling basis while prison interactions outside isolated areas continue as normal.
  
#Pison lockdown model
# The PPELockPrison folder contains the Prison isolation model. Use the individual.o, main.o, and individual.c files here to run the baseline model as defined below:
# Symptomatic inmates are tested using PCR with results by the next day. Prisons with symptomatic inmates with confirmed COVID-19 infection are locked down until no one is actively infected with COVID-19; N95 masks are used by staff interacting with isolated inmates; standard PPE masks are used by isolated inmates; prisons are locked on a rolling basis; inmates within the centre are free to move within the prison but there are no movements to or from other prisons (which operate as normal).
# Variations of this model implementing a 1-week, 3-week, 6-week, and 9-week delay in initiating prison lockdown after first detection of COVD-19 case can be found inside the PPELockPrison folder.

#Immunisation model
# The Immunisation folder contains the Immunisation model. Use the individual.o, main.o, and individual.c files here to run the baseline model as defined below:
# Assumes standard PPE is in place; 80% of the inmate population are assumed to have had a single dose and 50% to have had double dose immunisation; the same vaccination rate is applied for correctional and healthcare staff. We assumed the use of Pfizer mRNA vaccine. After the first dose, we applied a 46% reduction in onward COVID-19 transmission a 30% reduction in COVID-19 infection. After two doses, we applied a 65% reduction in onward COVID-19 transmission a 79% reduction in COVID-19 infection.
# Variations of this model implementing a high coverage or low coverage for inmates and staff are found inside the Immunisation folder.


  
  



  

 

  

  
