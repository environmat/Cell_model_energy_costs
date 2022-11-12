import cell_model_v004 as Cellmodel
# Active materials
NMC_cathode = Cellmodel.Activematerial_cathode("NMC811" , 195, 3.86)
LFP_cathode = Cellmodel.Activematerial_cathode("LFP", 160, 3.3)
Gr = Cellmodel.Activematerial_anode("Gr", 344 , 0.17)

# Currentcollectors 
Al = Cellmodel.CurrentCollector_cathode("Al", 14e-4 , 2.76)
Cu = Cellmodel.CurrentCollector_anode("Cu", 8e-4 , 8.96)

# Electrolyte 
LP40_standard = Cellmodel.Electrolyte("LiPF6" , "EC:DMC", 1.1)

# Separator
Separator_standard = Cellmodel.Separator("PP+Al" , 12e-4 , 0.44 , 1.18)

# Electrods 
# positive = Cellmodel.Electrodecomposition_cathode_opt1( NMC_cathode.name , 3.3 , 0.95 , 3.4, NMC_cathode) # NMC = 3.4 / LFP = 2.5
positive = Cellmodel.Electrodecomposition_cathode_opt1( LFP_cathode.name , 3.3 , 0.95 , 2.5, LFP_cathode) # NMC = 3.4 / LFP = 2.5
negative = Cellmodel.Electrodecomposition_anode_opt1( Gr.name , 3.3*1.1 , 0.965 , 1.6, Gr)

### Total cells
# test_cell = Cellmodel.Cylindrical("NMC-Gr_Cyl", positive, negative, Separator_standard , LP40_standard, 
                              # 1.7, 2.1, 7.0, 0.0165, 7.9, 0.25, 0.6, 4, 0.94, NMC_cathode, Gr, Al , Cu)
test_cell = Cellmodel.Cylindrical("LFP-Gr_Cyl", positive, negative, Separator_standard , LP40_standard, 
                              1.7, 2.1, 7.0, 0.0165, 7.9, 0.25, 0.6, 4, 0.94, NMC_cathode, Gr, Al , Cu)

 
 
# Masses
cellmass, Housing_mass = Cellmodel.getMass_cylindrical_total(positive, negative, Al, Cu, Separator_standard, test_cell, LP40_standard)

## Cathode Masses
# Cylindrical
cathode_mass = Cellmodel.getMass_cathode(positive) * test_cell.jr_area * 2


## Seperator and Elyte
Separator_mass = Cellmodel.getMass_separator(Separator_standard)
Elyte_mass = Cellmodel.getMass_electrolyte(LP40_standard, test_cell)

### Materials
## Cylindrical
Ni, Co, Mn, Li, Al = Cellmodel.getMass_elements_cyl_c(positive, Al, NMC_cathode, test_cell)
Gr, Cu, Binder_a = Cellmodel.getMass_elements_cyl_a(negative, Cu, test_cell)
Li_elyte = Cellmodel.getMass_elements_cyl_elyte(LP40_standard, Separator_standard, test_cell)

## Specific Lithium costs and processing
Pro_fac = 5.3 # CellEst Processing Factor Li_CO3
Li_CO3 = 53500 # $/mt https://www.spglobal.com/platts/en/market-insights/latest-news/metals/020322-chinese-lithium-carbonate-hydroxide-price-spread-at-record-high
LiOH = 47000 # $/mt

Li_CO3_kg = Li_CO3/1000
LiOH_kg = LiOH/1000

Li_LFP = Li_CO3_kg * Pro_fac
Li_NMC = LiOH_kg * Pro_fac

#### Material Costs
Ni_raw = 16
Co_raw = 51
Mn_raw = 2
Li_raw = Li_LFP 
Al_raw = 2
Al_cc_raw = 6
Gr_raw = 12
Cu_raw = 9
Binder_raw = 10 
Elyte_raw = 15
Separator_raw = 80
Steel_raw = 2
Conductive_raw = 7

### Cylindrical 
# # Costs Cylindrical
Ni_costs, Co_costs, Mn_costs, Li_costs, Al_costs, Al_cc_costs, Gr_costs, Cu_costs, \
    Binder_costs, Elyte_costs, Separator_costs, Housing_costs, Conductive_costs = \
    Cellmodel.getCosts_cyl(Ni_raw, Co_raw, Mn_raw, Li_raw, Al_raw, Al_cc_raw, Gr_raw, \
                           Cu_raw, Binder_raw, Elyte_raw, Separator_raw, Steel_raw, Conductive_raw)
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
print("Li content in Electrolyte [g]:", round((Li_elyte),2))
