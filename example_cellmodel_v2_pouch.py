import cell_model_v004 as Cellmodel
# Active materials
NMC_cathode = Cellmodel.Activematerial_cathode("NMC811" , 195, 3.86)
LFP_cathode = Cellmodel.Activematerial_cathode("LFP", 170, 3.3)
Gr = Cellmodel.Activematerial_anode("Gr", 344 , 0.17)

# Currentcollectors 
Al = Cellmodel.CurrentCollector_cathode("Al", 14e-4 , 2.76)
Cu = Cellmodel.CurrentCollector_anode("Cu", 8e-4 , 8.96)

# Electrolyte 
LP40_standard = Cellmodel.Electrolyte("LiPF6" , "EC:DMC", 1.1)

# Separator
Separator_standard = Cellmodel.Separator("PP+Al" , 12e-4 , 0.44 , 1.18)

# Electrods 
positive = Cellmodel.Electrodecomposition_cathode_opt1( NMC_cathode.name , 3.3 , 0.95 , 3.4, NMC_cathode) # NMC = 3.4 / LFP = 2.5
negative = Cellmodel.Electrodecomposition_anode_opt1( Gr.name , 3.3*1.1 , 0.965 , 1.6, Gr)

### Total cells
# NMC-Gr
test_cell = Cellmodel.Pouch("NMC-Gr" , positive, negative, Separator_standard , LP40_standard, 
                            1.7 , 30 , 10 , 35 , 113e-4, 1.8 , 1.0, 2 , 5 ,
                            15e-4, 2.76 , 8.96, 10, 0.93, NMC_cathode, Gr)
# LFP-Gr
# test_cell = Cellmodel.Pouch("LFP-Gr" , positive, negative, Separator_standard , LP40_standard, 
#                             1.7 , 9.5 , 6.1 , 7 , 113e-4, 1.8 , 1.0, 2 , 3 ,
#                             15e-4, 2.76 , 8.96, 4, 0.93, LFP_cathode, Gr)

 
 
# Masses
cellmass, pouch_packaging_mass = Cellmodel.getMass_Pouch_total(positive, negative, Al, Cu, Separator_standard, test_cell, LP40_standard)

 
## Cathode Masses
# Pouch
cathode_mass = Cellmodel.getMass_cathode(positive) * test_cell.jr_area * test_cell.nlayers

## Seperator and Elyte
Separator_mass = Cellmodel.getMass_separator(Separator_standard)
Elyte_mass = Cellmodel.getMass_electrolyte(LP40_standard, test_cell)

### Materials
## Pouch
Ni, Co, Mn, Li, Al = Cellmodel.getMass_elements_pouch_c(positive, Al, NMC_cathode, test_cell)
Gr, Cu, Binder_a = Cellmodel.getMass_elements_pouch_a(negative, Cu, test_cell) 


#### Material Costs
### Pouch
Ni_costs, Co_costs, Mn_costs, Li_costs, Al_costs, Al_cc_costs, Gr_costs, Cu_costs, \
    Binder_costs, Elyte_costs, Separator_costs, Housing_costs, Conductive_costs = \
    Cellmodel.getCosts_pouch(16, 51, 2, 45, 2, 6, 12, 9, 10, 15, 80, 2, 7)
CAM_Metal_costs = Ni_costs + Co_costs + Mn_costs + Li_costs + Al_costs# $ 
CAM_Metal_costs_kg = CAM_Metal_costs / (cathode_mass / 1000) # $ / kg 
Material_costs = Ni_costs + Co_costs + Mn_costs + Li_costs + Al_costs + Al_cc_costs + Gr_costs + Cu_costs + Binder_costs + \
                Elyte_costs + Separator_costs + Housing_costs + Conductive_costs


### All other costs in kWh
Separator_costs_kwh = Separator_costs / (test_cell.energy /1000)
Elyte_costs_kwh = Elyte_costs / (test_cell.energy /1000)
Al_cc_costs_kwh = Al_cc_costs / (test_cell.energy /1000)
Cu_costs_kwh = Cu_costs / (test_cell.energy /1000) 
Anode_costs_kwh = Gr_costs / (test_cell.energy /1000)
Binder_costs_kwh = Binder_costs / (test_cell.energy /1000)
Housing_costs_khw = Housing_costs / (test_cell.energy /1000)
Conductive_costs_kwh = Conductive_costs / (test_cell.energy /1000)

# Printing Info: 
print("Capacity [Ah]: ", round(test_cell.capacity,2))
print("Energy [Wh]: ", round(test_cell.energy,2))
print("Energy [Wh/kg]: ", round(test_cell.energy / (cellmass/1000),2))
print("CAM Metal costs [EUR/kWh]: ", round(CAM_Metal_costs / (test_cell.energy /1000),2))
print("AAM Metal costs [EUR/kWh]: ", round(Gr_costs / (test_cell.energy /1000),2))
print("Total Material costs [EUR/kWh]: ", round(Material_costs / (test_cell.energy /1000),2))
