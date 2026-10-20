rankine cycle simulator with irreversibilities in py (CoolProp)

the program prompts for the condenser and boiler pressures, the turbine inlet temperature and the isentropic efficencies of the pump and turbines.

model assumptions:
- state 1: saturated liquid at the condenser pressure
- pump: incompressible liquid
- boiler and condenser: isobaric processes, with no pressure drop
- turbine and pump: adiabatic, with isentropic efficiency
- energy changes are neglected

libraries:
- pip install CoolProp numpy matplotlib

execuction test (based on literature):

  - input parameters:
    condenser pressure: 0.75 bar
    boiler pressure: 30 bar
    turbine inlet temperature: 350 C
    efficencies: 1

  - terminal output summary:
states      T (°C)   P (bar)   h (kJ/kg)   s (kJ/kg.K)
1            91.76      0.75      384.44        1.2132
2            91.94     30.00      387.47        1.2132
3           350.00     30.00     3116.06        6.7449
4            91.76      0.75     2403.02        6.7449

specific pump work: 3.03 kJ/kg
specific turbine work: 713.04 kJ/kg
net specific work: 710.01 kJ/kg
heat added in boiler: 2728.59 kJ/kg
heat rejected in condenser: 2018.58 kJ/kg
overall thermal efficiency: 26.02 %
carnot efficiency (limit): 41.44 % 
steam quality at turbine exit: 88.61 %

  - literature validation:
    
  


