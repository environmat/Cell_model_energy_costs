## Source: 
# https://github.com/mjlacey/cellmodels/blob/master/cellmodel.jl

import scipy.integrate as integrate
import math

## Class definitions
## ============================================================================
class Activematerial_cathode:
    def __init__(self, name, spec_cap, Avg_E):
        self.name = name                    #           NMC811 / NMC622 / NCA / LFP
        self.spec_cap = spec_cap            # mAh g-1   195 / 181 / 210 / 170 
        self.Avg_E = Avg_E                  # V         3.86 / 3.86 / 3.86 / 3.3 
        
class Activematerial_anode:
    def __init__(self, name, spec_cap, Avg_E):
        self.name = name                # Gr / GrSi3_5
        self.spec_cap = spec_cap # mAh g-1   # 381 / 445 
        self.Avg_E = Avg_E # V            # 0.17 / 0.18

class CurrentCollector_cathode():
    def __init__(self, name_cc, thickness_cc_c, density_cc_c):
        self.name_cc = name_cc                  # Al / Cu
        self.thickness_cc_c = thickness_cc_c    # cm    # 8e-4 / 14e-4        
        self.density_cc_c = density_cc_c        # g cm-3    # 2.76 / 8.96
        self.mass_cc_c = self.thickness_cc_c * self.density_cc_c # g / cm2

class CurrentCollector_anode():
    def __init__(self, name_cc, thickness_cc_a, density_cc_a):
        self.name_cc = name_cc               # Al / Cu
        self.thickness_cc_a = thickness_cc_a # m    # 8e-6 / 14e-6        
        self.density_cc_a = density_cc_a # g cm-3    # 2.76 / 8.96
        self.mass_cc_a = self.thickness_cc_a * self.density_cc_a
# class Electrodecomposition_cathode_opt1:  ##### Activematerial_cathode.spec_cap does not work
#     def __init__(self, active_material, 
#                  areal_cap, active_frac, density_el):
#         self.active_material = active_material
#         self.areal_cap = areal_cap # mAh cm-2
#         self.active_frac = active_frac # Percent 
#         self.density_el = density_el # g cm-3
#     @property
#     def thickness_el(self, Activematerial_cathode): 
#         self.thickness_el_c = self.areal_cap / (self.active_frac * Activematerial_cathode.spec_cap * self.density_el)
#         return self.thickness_el_c
#     @property
#     def active_load(self, Activematerial_cathode):
#         self.active_load = self.areal_cap / Activematerial_cathode.spec_cap
#         return self.active_load
    
####
class Electrodecomposition_cathode_opt1():  
    def __init__(self, active_material, 
                 areal_cap, active_frac, density_el, Activematerial_cathode):
        self.active_material = active_material
        self.areal_cap = areal_cap # mAh cm-2
        self.active_frac = active_frac # Percent 
        self.density_el = density_el # g cm-3
        self.thickness_el_c = self.areal_cap / (self.active_frac * Activematerial_cathode.spec_cap * self.density_el)
        self.active_load = self.areal_cap / Activematerial_cathode.spec_cap
####

class Electrodecomposition_cathode_opt2: 
    def __init__(self, active_material, thickness_el_c, 
                 active_load, active_frac):
        self.active_material = active_material
        self.thickness_el_c = thickness_el_c # m
        self.active_load = active_load #mg cm-2
        self.active_frac = active_frac # Percent 
    @property     
    def areal_cap(self, Activematerial_cathode):
        self._areal_cap = Activematerial_cathode.spec_cap * self.active_load
        return self._areal_cap
    @property
    def density_el(self):
        self._density_el = (self.active_load / self.active_frac) / self.thickness_el_c
        return self._density_el

# class Electrodecomposition_anode_opt1: 
#     def __init__(self, active_material, 
#                  areal_cap, active_frac, density_el):
#         self.active_material = active_material
#         self.areal_cap = areal_cap # mAh cm-2
#         self.active_fraq = active_frac # Percent 
#         self.density_el = density_el # g cm-3
#     @property
#     def thickness_el(self, Activematerial_anode):
#         self._thickness_el_a = self.areal_cap / (self.active_frac * Activematerial_anode.spec_cap * self.density_el)
#         return self._thickness_el_a
#     @property
#     def active_load(self, Activematerial_anode):    
#         self._active_load = self.areal_cap / Activematerial_anode.spec_cap
#         return self._active_load

###
class Electrodecomposition_anode_opt1(): 
    def __init__(self, active_material, 
                 areal_cap, active_frac, density_el, Activematerial_anode):
        self.active_material = active_material
        self.areal_cap = areal_cap # mAh cm-2
        self.active_frac = active_frac # Percent 
        self.density_el = density_el # g cm-3
        self.thickness_el_a = self.areal_cap / (self.active_frac * Activematerial_anode.spec_cap * self.density_el)  
        self.active_load = self.areal_cap / Activematerial_anode.spec_cap
###

class Electrodecomposition_anode_opt2: 
    def __init__(self, active_material, thickness_el_a, 
                 active_load, active_frac):
        self.active_material = active_material
        self.thickness_el_a = thickness_el_a # m
        self.active_load = active_load #mg cm-2
        self.active_frac = active_frac # Percent 
    @property     
    def areal_cap(self, Activematerial_anode):
        self._areal_cap = Activematerial_anode.spec_cap * self.active_load
        return self._areal_cap
    @property
    def density_el(self):
        self._density_el = (self.active_load / self.active_frac) / self.thickness_el_a
        return self._density_el 

class Electrode:
    def __init__(self, Composite, CC):
        self.Composite = Composite
        self.CC = CC
            
class Separator:
    def __init__(self, name_sep, thickness_sep, porositiy_sep, density_sep):
        self.name_sep = name_sep
        self.thickness_sep = thickness_sep # 20 microns
        self.porositiy_sep = porositiy_sep
        self.density_sep = density_sep 

class Electrolyte_liquid:
    def __init__(self, salt, solvent, concentration, molarmass_salt):
        self.salt = salt
        self.solvent = solvent
        self.concentration = concentration
        self.molarmass_salt = molarmass_salt
    @property
    def saltmassfrac(self):
        self._saltmassfrac = 0.1222 * self.concentration
        return self._saltmassfrac
    @property
    def density_elyte(self):
        self._density_elyte = 0.7641*self.saltmassfrac+1.1299 # 1.1299 == density EC:DEC g cm-3
        return self._density_elyte

class Electrolyte_solid:
    def __init__(self, name_elyte, molarmass_el, molarmass_Li, index_Li, pos_electrode_fraction, density_elyte):
        self.name = name_elyte
        self.molarmass_el = molarmass_el # molarmass in g/mol
        self.molarmass_Li = molarmass_Li # molarmass in g/mol
        self.index_Li = index_Li # index in SE formula of Li, x in i.e. Li_x La_y Zr_z O_w
        self.pos_electrode_fraction = pos_electrode_fraction # share of pos electrode
        self.density_elyte = density_elyte # density in g cm-3

## Pouch ##
## ============================================================================
class Pouch:
    def __init__(self, name_pouch, positive, negative, separator, electrolyte, ecap_ratio, height, width, nlayers,
                  pouchthickness, pouchdensity, pouchclearance, terminalheight, terminalwidth, terminalthickness,
                  terminaldensitypos, terminaldensityneg, extramass, llifactor, 
                  Activematerial_cathode, Activematerial_anode):
        self.name_pouch = name_pouch
        self.positive = positive
        self.negative = negative
        self.separator = separator
        self.electrolyte = electrolyte
        self.ecap_ratio = ecap_ratio
        self.height = height
        self.width = width
        self.nlayers = nlayers
        self.pouchthickness = pouchthickness
        self.pouchdensity = pouchdensity
        self.pouchclearance = pouchclearance
        self.terminalheight = terminalheight
        self.terminalwidth = terminalwidth
        self.terminalthickness = terminalthickness
        self.terminaldensitypos = terminaldensitypos
        self.terminaldensityneg = terminaldensityneg
        self.extramass = extramass
        self.llifactor = llifactor
        self.jr_area = self.nlayers * self.width * self.height
        self.capacity = min(positive.areal_cap, negative.areal_cap) * self.llifactor * self.jr_area * 2 / 1000
        self.energy = (Activematerial_cathode.Avg_E - Activematerial_anode.Avg_E) * self.capacity


## Cylindrical ##
## ============================================================================
class Cylindrical:
    def __init__(self, name_cylindrical, positive, negative, separator, electrolyte, ecap_ratio, diameter, 
                 height, canthickness, candensity, voiddiameter, headspace, extramass, llifactor, 
                 Activematerial_cathode, Activematerial_anode, CurrentCollector_cathode,
                 CurrentCollector_anode):
        self.name_cylindrical = name_cylindrical
        self.positive = positive
        self.negative = negative
        self.separator = separator
        self.electrolyte = electrolyte
        self.ecap_ratio = ecap_ratio
        self.diameter = diameter
        self.height = height
        self.canthickness = canthickness
        self.candensity = candensity
        self.voiddiameter = voiddiameter
        self.headspace = headspace
        self.extramass = extramass
        self.llifactor = llifactor
        self.stackthickness = (2 * positive.thickness_el_c) + (2 * negative.thickness_el_a) + (2 * separator.thickness_sep) + \
                              CurrentCollector_cathode.thickness_cc_c + CurrentCollector_anode.thickness_cc_a
        self.turns = (self.diameter/2 - self.voiddiameter/2 - self.canthickness - self.stackthickness ) / self.stackthickness
                        # ((self.diameter - (2 * self.canthickness) - self.stackthickness - 
                        # self.voiddiameter) / 2) / self.stackthickness
                     
        def integrand(x):
            return math.sqrt(((self.voiddiameter / 2) + self.stackthickness * x / (2*math.pi))**2 + \
                   (self.stackthickness / (2*math.pi))**2)
        self.length = integrate.quad(integrand, 0, self.turns * 2 * math.pi)[0]
        self.jr_area = self.length * (self.height - self.headspace - 2 * self.canthickness)
        self.capacity = min(positive.areal_cap, negative.areal_cap) * self.llifactor * self.jr_area * 2 / 1000
        self.energy = (Activematerial_cathode.Avg_E - Activematerial_anode.Avg_E) * self.capacity


## Prismatic "Cinnamonroll" ##
## ============================================================================
class Prismatic:
    def __init__(self, name_prismatic, positive, negative, separator, electrolyte, ecap_ratio,  
                 height, width, depth, canthickness, candensity, termclearance, nrolls, extramass, llifactor, 
                 Activematerial_cathode, Activematerial_anode, CurrentCollector_cathode,
                 CurrentCollector_anode):
        self.name_prismatic = name_prismatic
        self.positive = positive
        self.negative = negative
        self.separator = separator
        self.electrolyte = electrolyte
        self.ecap_ratio = ecap_ratio
        self.height = height
        self.width = width
        self.depth = depth
        self.canthickness = canthickness
        self.candensity = candensity
        self.termclearance = termclearance
        self.nrolls = nrolls
        self.extramass = extramass
        self.llifactor = llifactor
        self.stackthickness = (2 * positive.thickness_el_c) + (2 * negative.thickness_el_a) + (2 * separator.thickness_sep) + \
                              CurrentCollector_cathode.thickness_cc_c + CurrentCollector_anode.thickness_cc_a
        self.rollthickness = (self.depth - (2 * self.canthickness) - 120e-4) / self.nrolls
        self.turns = (( self.rollthickness - self.stackthickness - 2 * separator.thickness_sep) / 2 ) / self.stackthickness
        def integrand(x):
            return math.sqrt((self.stackthickness * x / (2*math.pi))**2 + (self.stackthickness / (2*math.pi))**2)
        self.length = integrate.quad(integrand, 0, self.turns * 2 * math.pi)[0]
        self.jr_area = (self.length + ((self.width - self.rollthickness) * 2 * math.floor(self.turns))) * \
                       (self.height - self.termclearance) * self.nrolls
        self.capacity = min(positive.areal_cap, negative.areal_cap) * self.llifactor * self.jr_area * 2 / 1000
        self.energy = (Activematerial_cathode.Avg_E - Activematerial_anode.Avg_E) * self.capacity

        

## Masses ##
### ===========================================================================
def getMass_cathode(Electrodecomposition_cathode_opt1):
    mass_cathode = Electrodecomposition_cathode_opt1.active_load / Electrodecomposition_cathode_opt1.active_frac     
    return mass_cathode

def getMass_anode(Electrodecomposition_anode_opt1):
    mass_anode = Electrodecomposition_anode_opt1.active_load / Electrodecomposition_anode_opt1.active_frac    
    return mass_anode

def getMass_separator(Separator):
    getMass_separator.mass_sep = Separator.thickness_sep * Separator.density_sep # g / cm2
    return getMass_separator.mass_sep

def getMass_cc_c(CurrentCollector_cathode):
    mass_cc_c = CurrentCollector_cathode.mass_cc_c
    return mass_cc_c

def getMass_cc_a(CurrentCollector_anode):
    mass_cc_a = CurrentCollector_anode.mass_cc_a
    return mass_cc_a

def getMass_electrode_c(Electrodecomposition_cathode_opt1, CurrentCollector_cathode): 
    mass_cathode = Electrodecomposition_cathode_opt1.active_load / Electrodecomposition_cathode_opt1.active_frac
    mass_cc_c = CurrentCollector_cathode.mass_cc_c
    mass_elec_c = 2 * mass_cathode + mass_cc_c
    return mass_elec_c

def getMass_electrode_a(Electrodecomposition_anode_opt1, CurrentCollector_anode): 
    mass_anode = Electrodecomposition_anode_opt1.active_load / Electrodecomposition_anode_opt1.active_frac
    mass_cc_a = CurrentCollector_anode.mass_cc_a
    mass_elec_a = 2 * mass_anode + mass_cc_a
    return mass_elec_a

def getMass_stack(Electrodecomposition_cathode_opt1, Electrodecomposition_anode_opt1, \
                  CurrentCollector_cathode, CurrentCollector_anode, Separator):
    #
    mass_cathode = Electrodecomposition_cathode_opt1.active_load / Electrodecomposition_cathode_opt1.active_frac
    mass_cc_c = CurrentCollector_cathode.mass_cc_c
    mass_elec_c = 2 * mass_cathode + mass_cc_c
    #
    mass_anode = Electrodecomposition_anode_opt1.active_load / Electrodecomposition_anode_opt1.active_frac
    mass_cc_a = CurrentCollector_anode.mass_cc_a
    mass_elec_a = 2 * mass_anode + mass_cc_a
    #
    mass_sep = Separator.thickness_sep * Separator.density_sep
    #     
    mass_stack = mass_elec_c + mass_elec_a + mass_sep 
    return mass_stack

def getMass_jr_cyl(Electrodecomposition_cathode_opt1, Electrodecomposition_anode_opt1, \
                  CurrentCollector_cathode, CurrentCollector_anode, Separator, Cylindrical):
    #
    mass_cathode = Electrodecomposition_cathode_opt1.active_load / Electrodecomposition_cathode_opt1.active_frac
    mass_cc_c = CurrentCollector_cathode.mass_cc_c
    mass_elec_c = 2 * mass_cathode + mass_cc_c
    #
    mass_anode = Electrodecomposition_anode_opt1.active_load / Electrodecomposition_anode_opt1.active_frac
    mass_cc_a = CurrentCollector_anode.mass_cc_a
    mass_elec_a = 2 * mass_anode + mass_cc_a
    #
    mass_sep = Separator.thickness_sep * Separator.density_sep
    #     
    mass_stack = mass_elec_c + mass_elec_a + mass_sep 
    #
    mass_jr = mass_stack * Cylindrical.jr_area
    return mass_jr

def getMass_electrolyte_liquid(Electrolyte_liquid, Cylindrical):
    getMass_electrolyte_liquid.mass_elyte = Electrolyte_liquid.density_elyte * Cylindrical.ecap_ratio * Cylindrical.capacity
    getMass_electrolyte_liquid.mass_elyte_salt = Electrolyte_liquid.molarmass_salt * (Cylindrical.ecap_ratio/1000) * Cylindrical.capacity * Electrolyte_liquid.concentration
                                        # (g/mol)*(mL/Ah)*(Ah)*(mol/L)
    return getMass_electrolyte_liquid.mass_elyte, getMass_electrolyte_liquid.mass_elyte_salt

def getMass_cylindrical_housing(Cylindrical): 
    vol_cylindrical_housing = (math.pi * (Cylindrical.diameter / 2)**2 - math.pi * ((Cylindrical.diameter - 
                               Cylindrical.canthickness)/2)**2) * Cylindrical.height + \
                               (math.pi * (Cylindrical.diameter / 2)**2) * Cylindrical.canthickness
    mass_cylindrical_housing = vol_cylindrical_housing * Cylindrical.candensity + Cylindrical.extramass
    return mass_cylindrical_housing

#######################
# Sum up cylindrical masses 
def getMass_cylindrical_total(Electrodecomposition_cathode_opt1, Electrodecomposition_anode_opt1, \
                             CurrentCollector_cathode, CurrentCollector_anode, Separator, Cylindrical, \
                             Electrolyte):
    #
    mass_cathode = Electrodecomposition_cathode_opt1.active_load / Electrodecomposition_cathode_opt1.active_frac
    mass_cc_c = CurrentCollector_cathode.mass_cc_c
    mass_elec_c = 2 * mass_cathode + mass_cc_c
    #
    mass_anode = Electrodecomposition_anode_opt1.active_load / Electrodecomposition_anode_opt1.active_frac
    mass_cc_a = CurrentCollector_anode.mass_cc_a
    mass_elec_a = 2 * mass_anode + mass_cc_a
    #
    mass_sep = Separator.thickness_sep * Separator.density_sep
    #     
    mass_stack = mass_elec_c + mass_elec_a + mass_sep 
    #
    mass_jr = mass_stack * Cylindrical.jr_area
    #
    mass_elyte = Electrolyte.density_elyte * Cylindrical.ecap_ratio * Cylindrical.capacity
    #
    vol_cylindrical_housing = (math.pi * (Cylindrical.diameter / 2)**2 - math.pi * ((Cylindrical.diameter - 
                               Cylindrical.canthickness)/2)**2) * Cylindrical.height + \
                               (math.pi * (Cylindrical.diameter / 2)**2) * Cylindrical.canthickness
    getMass_cylindrical_total.mass_cylindrical_housing = vol_cylindrical_housing * Cylindrical.candensity + Cylindrical.extramass
    #
    getMass_cylindrical_total.mass_total = mass_jr + mass_elyte + getMass_cylindrical_total.mass_cylindrical_housing + \
                                            Cylindrical.extramass
    return getMass_cylindrical_total.mass_total , getMass_cylindrical_total.mass_cylindrical_housing

    
def getMass_pouch(Pouch):
     pouchmass = Pouch.height * Pouch.width * Pouch.pouchthickness * Pouch.pouchdensity
     
     terminalmass = (Pouch.terminalheight * Pouch.terminalwidth * Pouch.terminalthickness) * \
                    (Pouch.terminaldensitypos + Pouch.terminaldensityneg) 
     return pouchmass 

#######################
# Sum up pouch masses  
def getMass_Pouch_total(Electrodecomposition_cathode_opt1, Electrodecomposition_anode_opt1, \
                             CurrentCollector_cathode, CurrentCollector_anode, Separator, Pouch, \
                             Electrolyte):
    #
    mass_cathode = Electrodecomposition_cathode_opt1.active_load / Electrodecomposition_cathode_opt1.active_frac
    mass_cc_c = CurrentCollector_cathode.mass_cc_c
    mass_elec_c = 2 * mass_cathode + mass_cc_c
    #
    mass_anode = Electrodecomposition_anode_opt1.active_load / Electrodecomposition_anode_opt1.active_frac
    mass_cc_a = CurrentCollector_anode.mass_cc_a
    mass_elec_a = 2 * mass_anode + mass_cc_a
    #
    mass_sep = Separator.thickness_sep * Separator.density_sep
    #     
    mass_stack = mass_elec_c + mass_elec_a + mass_sep 
    #
    mass_stack_nlayers = mass_stack * Pouch.jr_area
    #
    mass_elyte = Electrolyte.density_elyte * Pouch.ecap_ratio * Pouch.capacity
    #
    pouchmass = (Pouch.height + Pouch.pouchclearance) * (Pouch.width + Pouch.pouchclearance) * Pouch.pouchthickness * Pouch.pouchdensity * 2

    terminalmass = (Pouch.terminalheight + Pouch.pouchclearance) * (Pouch.terminalwidth + Pouch.pouchclearance) * Pouch.terminalthickness * Pouch.terminaldensitypos + \
        (Pouch.terminalheight + Pouch.pouchclearance) * (Pouch.terminalwidth) * Pouch.terminalthickness * Pouch.terminaldensityneg

    #
    getMass_Pouch_total.mass_packaging_pouch = pouchmass + terminalmass + Pouch.extramass
    #
    getMass_Pouch_total.mass_total = mass_stack_nlayers + mass_elyte + getMass_Pouch_total.mass_packaging_pouch
    #
    return getMass_Pouch_total.mass_total , getMass_Pouch_total.mass_packaging_pouch

#######################
# Sum up prismatic masses  
def getMass_prismatic_total(Electrodecomposition_cathode_opt1, Electrodecomposition_anode_opt1, \
                             CurrentCollector_cathode, CurrentCollector_anode, Separator, Prismatic, \
                             Electrolyte):
    #
    mass_cathode = Electrodecomposition_cathode_opt1.active_load / Electrodecomposition_cathode_opt1.active_frac
    mass_cc_c = CurrentCollector_cathode.mass_cc_c
    mass_elec_c = 2 * mass_cathode + mass_cc_c
    #
    mass_anode = Electrodecomposition_anode_opt1.active_load / Electrodecomposition_anode_opt1.active_frac
    mass_cc_a = CurrentCollector_anode.mass_cc_a
    mass_elec_a = 2 * mass_anode + mass_cc_a
    #
    mass_sep = Separator.thickness_sep * Separator.density_sep
    #     
    mass_stack = mass_elec_c + mass_elec_a + mass_sep 
    #
    mass_stack_nlayers = mass_stack * Prismatic.jr_area
    #
    mass_elyte = Electrolyte.density_elyte * Prismatic.ecap_ratio * Prismatic.capacity
    #
    prismatic_housing = ((Prismatic.width * Prismatic.height * Prismatic.depth) - \
        ((Prismatic.width - Prismatic.canthickness) * (Prismatic.height - Prismatic.canthickness) * \
         (Prismatic.depth - Prismatic.canthickness))) * Prismatic.candensity 
    #
    getMass_prismatic_total.mass_packaging_prismatic = prismatic_housing + Prismatic.extramass
    #
    getMass_prismatic_total.mass_total = mass_stack_nlayers + mass_elyte + getMass_prismatic_total.mass_packaging_prismatic
    #
    return getMass_prismatic_total.mass_total , getMass_prismatic_total.mass_packaging_prismatic

############################
# Mass for cylindrical cell     
def getMass_elements_cyl_c(Electrodecomposition_cathode_opt1, CurrentCollector_cathode, Activematerial_cathode, Cylindrical):
    #
    mass_cathode = Electrodecomposition_cathode_opt1.active_load / Electrodecomposition_cathode_opt1.active_frac
    mass_cc_c = CurrentCollector_cathode.mass_cc_c
    # NMC811
    if Activematerial_cathode.name == "NMC811": # energies-434768-supplementary.xlsx
        getMass_elements_cyl_c.Ni_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.48 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_cyl_c.Co_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.06 * 0.9
        getMass_elements_cyl_c.Mn_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.06 * 0.9
        getMass_elements_cyl_c.Li_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.07 * 0.9
        getMass_elements_cyl_c.Al_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.00 * 0.9
        getMass_elements_cyl_c.cc_c_mass = mass_cc_c * Cylindrical.jr_area
        getMass_elements_cyl_c.binder_mass_c = (2 * mass_cathode * Cylindrical.jr_area) * 0.05 # 5% of total share
        getMass_elements_cyl_c.conductive_material = (2 * mass_cathode * Cylindrical.jr_area) * 0.05 # 5% of total share
        getMass_elements_cyl_c.Separator_mass = 2 * getMass_separator.mass_sep * Cylindrical.jr_area
    # NMC622
    elif Activematerial_cathode.name == "NMC622":
        getMass_elements_cyl_c.Ni_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.36 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_cyl_c.Co_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.12 * 0.9
        getMass_elements_cyl_c.Mn_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.12 * 0.9
        getMass_elements_cyl_c.Li_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.07 * 0.9
        getMass_elements_cyl_c.Al_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.00 * 0.9
        getMass_elements_cyl_c.cc_c_mass = mass_cc_c * Cylindrical.jr_area
        getMass_elements_cyl_c.binder_mass_c = (2 * mass_cathode * Cylindrical.jr_area) * 0.05 # 5% of total share
        getMass_elements_cyl_c.conductive_material = (2 * mass_cathode * Cylindrical.jr_area) * 0.05 # 5% of total share
        getMass_elements_cyl_c.Separator_mass = 2 * getMass_separator.mass_sep * Cylindrical.jr_area
    # NMC333
    elif Activematerial_cathode.name == "NMC333":
        getMass_elements_cyl_c.Ni_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.20 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_cyl_c.Co_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.20 * 0.9
        getMass_elements_cyl_c.Mn_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.20 * 0.9
        getMass_elements_cyl_c.Li_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.07 * 0.9
        getMass_elements_cyl_c.Al_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.00 * 0.9
        getMass_elements_cyl_c.cc_c_mass = mass_cc_c * Cylindrical.jr_area
        getMass_elements_cyl_c.binder_mass_c = (2 * mass_cathode * Cylindrical.jr_area) * 0.05 # 5% of total share
        getMass_elements_cyl_c.conductive_material = (2 * mass_cathode * Cylindrical.jr_area) * 0.05 # 5% of total share
        getMass_elements_cyl_c.Separator_mass = 2 * getMass_separator.mass_sep * Cylindrical.jr_area
    # NCA
    elif Activematerial_cathode.name == "NCA":
        getMass_elements_cyl_c.Ni_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.47 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_cyl_c.Co_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.09 * 0.9
        getMass_elements_cyl_c.Mn_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.00 * 0.9
        getMass_elements_cyl_c.Li_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.07 * 0.9
        getMass_elements_cyl_c.Al_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.01 * 0.9
        getMass_elements_cyl_c.cc_c_mass = mass_cc_c * Cylindrical.jr_area
        getMass_elements_cyl_c.binder_mass_c = (2 * mass_cathode * Cylindrical.jr_area) * 0.05 # 5% of total share
        getMass_elements_cyl_c.conductive_material = (2 * mass_cathode * Cylindrical.jr_area) * 0.05 # 5% of total share
        getMass_elements_cyl_c.Separator_mass = 2 * getMass_separator.mass_sep * Cylindrical.jr_area
    # LFP
    elif Activematerial_cathode.name == "LFP":
        getMass_elements_cyl_c.Ni_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.00 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_cyl_c.Co_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.00 * 0.9
        getMass_elements_cyl_c.Mn_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.00 * 0.9
        getMass_elements_cyl_c.Li_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.04 * 0.9
        getMass_elements_cyl_c.Al_mass = (2 * mass_cathode * Cylindrical.jr_area) * 0.00 * 0.9
        getMass_elements_cyl_c.cc_c_mass = mass_cc_c * Cylindrical.jr_area
        getMass_elements_cyl_c.binder_mass_c = (2 * mass_cathode * Cylindrical.jr_area) * 0.05 # 5% of total share
        getMass_elements_cyl_c.conductive_material = (2 * mass_cathode * Cylindrical.jr_area) * 0.05 # 5% of total share
        getMass_elements_cyl_c.Separator_mass = 2 * getMass_separator.mass_sep * Cylindrical.jr_area
    else:
        1
    return getMass_elements_cyl_c.Ni_mass, getMass_elements_cyl_c.Co_mass, \
           getMass_elements_cyl_c.Mn_mass, getMass_elements_cyl_c.Li_mass, \
           getMass_elements_cyl_c.Al_mass

# individual cell chemistries
def getMass_elements_cyl_c_spec(Electrodecomposition_cathode_opt1, CurrentCollector_cathode, Activematerial_cathode, Cylindrical,\
                                Ni_content, Co_content, Mn_content, Li_content, Al_content, Active_material_share, Binder_share, Conductive_material_share):
    #
    mass_cathode = Electrodecomposition_cathode_opt1.active_load / Electrodecomposition_cathode_opt1.active_frac
    mass_cc_c = CurrentCollector_cathode.mass_cc_c
    #
    getMass_elements_cyl_c_spec.Ni_mass = (2 * mass_cathode * Cylindrical.jr_area) * Ni_content * Active_material_share # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
    getMass_elements_cyl_c_spec.Co_mass = (2 * mass_cathode * Cylindrical.jr_area) * Co_content * Active_material_share
    getMass_elements_cyl_c_spec.Mn_mass = (2 * mass_cathode * Cylindrical.jr_area) * Mn_content * Active_material_share
    getMass_elements_cyl_c_spec.Li_mass = (2 * mass_cathode * Cylindrical.jr_area) * Li_content * Active_material_share
    getMass_elements_cyl_c_spec.Al_mass = (2 * mass_cathode * Cylindrical.jr_area) * Al_content * Active_material_share
    getMass_elements_cyl_c_spec.cc_c_mass = mass_cc_c * Cylindrical.jr_area
    getMass_elements_cyl_c_spec.binder_mass_c = (2 * mass_cathode * Cylindrical.jr_area) * Binder_share # 5% of total share
    getMass_elements_cyl_c_spec.conductive_material = (2 * mass_cathode * Cylindrical.jr_area) * Conductive_material_share # 5% of total share
    getMass_elements_cyl_c_spec.Separator_mass = 2 * getMass_separator.mass_sep * Cylindrical.jr_area
    
    return getMass_elements_cyl_c_spec.Ni_mass, getMass_elements_cyl_c_spec.Co_mass, \
           getMass_elements_cyl_c_spec.Mn_mass, getMass_elements_cyl_c_spec.Li_mass, \
           getMass_elements_cyl_c_spec.Al_mass

def getMass_elements_cyl_a(Electrodecomposition_anode_opt1, CurrentCollector_anode, Cylindrical): 
    #
    mass_anode = Electrodecomposition_anode_opt1.active_load / Electrodecomposition_anode_opt1.active_frac
    mass_cc_a = CurrentCollector_anode.mass_cc_a
    # Graphite
    getMass_elements_cyl_a.Gr_mass = (2 * mass_anode * Cylindrical.jr_area) * Electrodecomposition_anode_opt1.active_frac
    getMass_elements_cyl_a.cc_a_mass = mass_cc_a * Cylindrical.jr_area
    getMass_elements_cyl_a.binder_mass_a = (2 * mass_anode * Cylindrical.jr_area) * 0.03 # 3% of total share
    return getMass_elements_cyl_a.Gr_mass, getMass_elements_cyl_a.cc_a_mass, \
           getMass_elements_cyl_a.binder_mass_a

def getMass_elements_cyl_elyte(Electrolyte_liquid, Separator, Cylindrical):
    #
    getMass_electrolyte_liquid.mass_elyte = Electrolyte_liquid.density_elyte * Cylindrical.ecap_ratio * Cylindrical.capacity
    if Electrolyte_liquid.salt == "LiPF6":
        Li_mass_content_LiPF6 = 0.0457 # 4.57 weight% of Li content in LiPF6 / concentration in mol L-1
        getMass_elements_cyl_elyte.Li_mass = getMass_electrolyte_liquid.mass_elyte_salt * Li_mass_content_LiPF6
    # elif Electrolyte.salt == "LLZO":
    #     getMass_separator.mass_sep = Separator.thickness_sep * Separator.density_sep  # g / cm2
    #     Li_mass_content_LLZO = 0.0579  # 5.79 weight% of Li content in LLZO
    #     getMass_elements_cyl_elyte.Li_mass = (Electrolyte.concentration * getMass_electrolyte.mass_elyte + getMass_separator.mass_sep) \
    #                                          * Li_mass_content_LLZO
    else:
        getMass_elements_cyl_elyte.Li_mass = 1
    return getMass_elements_cyl_elyte.Li_mass

def getMass_elements_elyte_solid(Electrolyte_solid, Electrodecomposition_cathode_opt1):
    #
    mass_cathode = Electrodecomposition_cathode_opt1.areal_cap / Electrodecomposition_cathode_opt1.active_frac
    #
    getMass_elements_elyte_solid.Li_mass_fraction_elyte = (Electrolyte_solid.index_Li * Electrolyte_solid.molarmass_Li) / Electrolyte_solid.molarmass_el
    getMass_elements_elyte_solid.Li_mass_abs = mass_cathode * Electrolyte_solid.pos_electrode_fraction * getMass_elements_elyte_solid.Li_mass_fraction_elyte

    # molarmass_el, molarmass_Li, index_Li, pos_electrode_fraction
    return getMass_elements_elyte_solid.Li_mass_abs

#######################
#######################
# Get costs for all components Cylindrical
def getCosts_cyl(Ni_raw , Co_raw, Mn_raw, Li_raw, Al_raw, Al_cc_raw, Gr_raw, Cu_raw, \
                 Binder_raw, Elyte_raw, Separator_raw, Steel_raw, Conductive_raw, Electrolyte_liquid):
    # Cathode
    Ni_costs = getMass_elements_cyl_c.Ni_mass/1000 * Ni_raw # gram / 1000 * $ / kg = $ 
    Co_costs = getMass_elements_cyl_c.Co_mass/1000 * Co_raw
    Mn_costs = getMass_elements_cyl_c.Mn_mass/1000 * Mn_raw
    Li_costs = getMass_elements_cyl_c.Li_mass/1000 * Li_raw
    Al_costs = getMass_elements_cyl_c.Al_mass/1000 * Al_raw
    Al_cc_costs = getMass_elements_cyl_c.cc_c_mass/1000 * Al_cc_raw
    # Anode
    Gr_costs = getMass_elements_cyl_a.Gr_mass/1000 * Gr_raw
    Cu_costs = getMass_elements_cyl_a.cc_a_mass/1000 * Cu_raw
    # Binder for both
    Binder_costs = (getMass_elements_cyl_c.binder_mass_c + getMass_elements_cyl_a.binder_mass_a)/1000 * Binder_raw
    # Electrolyte
    if Electrolyte_liquid.salt == "LiPF6" :
        Elyte_costs = getMass_electrolyte_liquid.mass_elyte/1000 * Elyte_raw
    else:
        Elyte_costs = getMass_elements_elyte_solid.Li_mass_abs/1000 * Li_raw
    # Separator
    Separator_costs = getMass_elements_cyl_c.Separator_mass/1000 * Separator_raw
    # Housing
    Housing_costs = getMass_cylindrical_total.mass_cylindrical_housing/1000 * Steel_raw
    # Conductive Material
    Conductive_material_costs = getMass_elements_cyl_c.conductive_material/1000 * Conductive_raw
    return Ni_costs, Co_costs, Mn_costs, Li_costs, Al_costs, Al_cc_costs, Gr_costs, Cu_costs, \
        Binder_costs, Elyte_costs, Separator_costs, Housing_costs, Conductive_material_costs
    
    
############################
# Mass for pouch cell     
def getMass_elements_pouch_c(Electrodecomposition_cathode_opt1, CurrentCollector_cathode, Activematerial_cathode, Pouch):
    #
    mass_cathode = Electrodecomposition_cathode_opt1.active_load / Electrodecomposition_cathode_opt1.active_frac
    mass_cc_c = CurrentCollector_cathode.mass_cc_c
    # NMC811
    if Activematerial_cathode.name == "NMC811": # energies-434768-supplementary.xlsx
        getMass_elements_pouch_c.Ni_mass = (2 * mass_cathode * Pouch.jr_area) * 0.48 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_pouch_c.Co_mass = (2 * mass_cathode * Pouch.jr_area) * 0.06 * 0.9
        getMass_elements_pouch_c.Mn_mass = (2 * mass_cathode * Pouch.jr_area) * 0.06 * 0.9
        getMass_elements_pouch_c.Li_mass = (2 * mass_cathode * Pouch.jr_area) * 0.07 * 0.9
        getMass_elements_pouch_c.Al_mass = (2 * mass_cathode * Pouch.jr_area) * 0.00 * 0.9
        getMass_elements_pouch_c.cc_c_mass = mass_cc_c * Pouch.jr_area
        getMass_elements_pouch_c.binder_mass_c = (2 * mass_cathode * Pouch.jr_area) * 0.05 # 5% of total share
        getMass_elements_pouch_c.conductive_material = (2 * mass_cathode * Pouch.jr_area) * 0.05 # 5% of total share
        getMass_elements_pouch_c.Separator_mass = 2 * getMass_separator.mass_sep * Pouch.jr_area
    # NMC622
    elif Activematerial_cathode.name == "NMC622":
        getMass_elements_pouch_c.Ni_mass = (2 * mass_cathode * Pouch.jr_area) * 0.36 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_pouch_c.Co_mass = (2 * mass_cathode * Pouch.jr_area) * 0.12 * 0.9
        getMass_elements_pouch_c.Mn_mass = (2 * mass_cathode * Pouch.jr_area) * 0.12 * 0.9
        getMass_elements_pouch_c.Li_mass = (2 * mass_cathode * Pouch.jr_area) * 0.07 * 0.9
        getMass_elements_pouch_c.Al_mass = (2 * mass_cathode * Pouch.jr_area) * 0.00 * 0.9
        getMass_elements_pouch_c.cc_c_mass = mass_cc_c * Pouch.jr_area
        getMass_elements_pouch_c.binder_mass_c = (2 * mass_cathode * Pouch.jr_area) * 0.05 # 5% of total share
        getMass_elements_pouch_c.conductive_material = (2 * mass_cathode * Pouch.jr_area) * 0.05 # 5% of total share
        getMass_elements_pouch_c.Separator_mass = 2 * getMass_separator.mass_sep * Pouch.jr_area
    # NMC333
    elif Activematerial_cathode.name == "NMC333":
        getMass_elements_pouch_c.Ni_mass = (2 * mass_cathode * Pouch.jr_area) * 0.20 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_pouch_c.Co_mass = (2 * mass_cathode * Pouch.jr_area) * 0.20 * 0.9
        getMass_elements_pouch_c.Mn_mass = (2 * mass_cathode * Pouch.jr_area) * 0.20 * 0.9
        getMass_elements_pouch_c.Li_mass = (2 * mass_cathode * Pouch.jr_area) * 0.07 * 0.9
        getMass_elements_pouch_c.Al_mass = (2 * mass_cathode * Pouch.jr_area) * 0.00 * 0.9
        getMass_elements_pouch_c.cc_c_mass = mass_cc_c * Pouch.jr_area
        getMass_elements_pouch_c.binder_mass_c = (2 * mass_cathode * Pouch.jr_area) * 0.05 # 5% of total share
        getMass_elements_pouch_c.conductive_material = (2 * mass_cathode * Pouch.jr_area) * 0.05 # 5% of total share
        getMass_elements_pouch_c.Separator_mass = 2 * getMass_separator.mass_sep * Pouch.jr_area
    # NCA
    elif Activematerial_cathode.name == "NCA":
        getMass_elements_pouch_c.Ni_mass = (2 * mass_cathode * Pouch.jr_area) * 0.47 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_pouch_c.Co_mass = (2 * mass_cathode * Pouch.jr_area) * 0.09 * 0.9
        getMass_elements_pouch_c.Mn_mass = (2 * mass_cathode * Pouch.jr_area) * 0.00 * 0.9
        getMass_elements_pouch_c.Li_mass = (2 * mass_cathode * Pouch.jr_area) * 0.07 * 0.9
        getMass_elements_pouch_c.Al_mass = (2 * mass_cathode * Pouch.jr_area) * 0.01 * 0.9
        getMass_elements_pouch_c.cc_c_mass = mass_cc_c * Pouch.jr_area
        getMass_elements_pouch_c.binder_mass_c = (2 * mass_cathode * Pouch.jr_area) * 0.05 # 5% of total share
        getMass_elements_pouch_c.conductive_material = (2 * mass_cathode * Pouch.jr_area) * 0.05 # 5% of total share
        getMass_elements_pouch_c.Separator_mass = 2 * getMass_separator.mass_sep * Pouch.jr_area
    # LFP
    elif Activematerial_cathode.name == "LFP":
        getMass_elements_pouch_c.Ni_mass = (2 * mass_cathode * Pouch.jr_area) * 0.00 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_pouch_c.Co_mass = (2 * mass_cathode * Pouch.jr_area) * 0.00 * 0.9
        getMass_elements_pouch_c.Mn_mass = (2 * mass_cathode * Pouch.jr_area) * 0.00 * 0.9
        getMass_elements_pouch_c.Li_mass = (2 * mass_cathode * Pouch.jr_area) * 0.04 * 0.9
        getMass_elements_pouch_c.Al_mass = (2 * mass_cathode * Pouch.jr_area) * 0.00 * 0.9
        getMass_elements_pouch_c.cc_c_mass = mass_cc_c * Pouch.jr_area
        getMass_elements_pouch_c.binder_mass_c = (2 * mass_cathode * Pouch.jr_area) * 0.05 # 5% of total share
        getMass_elements_pouch_c.conductive_material = (2 * mass_cathode * Pouch.jr_area) * 0.05 # 5% of total share
        getMass_elements_pouch_c.Separator_mass = 2 * getMass_separator.mass_sep * Pouch.jr_area
    else:
        1
    return getMass_elements_pouch_c.Ni_mass, getMass_elements_pouch_c.Co_mass, \
           getMass_elements_pouch_c.Mn_mass, getMass_elements_pouch_c.Li_mass, \
           getMass_elements_pouch_c.Al_mass

# individual cell chemistries
def getMass_elements_pouch_c_spec(Electrodecomposition_cathode_opt1, CurrentCollector_cathode, Activematerial_cathode, Pouch,\
                                Ni_content, Co_content, Mn_content, Li_content, Al_content, Active_material_share, Binder_share, Conductive_material_share):
    #
    mass_cathode = Electrodecomposition_cathode_opt1.active_load / Electrodecomposition_cathode_opt1.active_frac
    mass_cc_c = CurrentCollector_cathode.mass_cc_c
    #
    getMass_elements_pouch_c_spec.Ni_mass = (2 * mass_cathode * Pouch.jr_area) * Ni_content * Active_material_share # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
    getMass_elements_pouch_c_spec.Co_mass = (2 * mass_cathode * Pouch.jr_area) * Co_content * Active_material_share
    getMass_elements_pouch_c_spec.Mn_mass = (2 * mass_cathode * Pouch.jr_area) * Mn_content * Active_material_share
    getMass_elements_pouch_c_spec.Li_mass = (2 * mass_cathode * Pouch.jr_area) * Li_content * Active_material_share
    getMass_elements_pouch_c_spec.Al_mass = (2 * mass_cathode * Pouch.jr_area) * Al_content * Active_material_share
    getMass_elements_pouch_c_spec.cc_c_mass = mass_cc_c * Pouch.jr_area
    getMass_elements_pouch_c_spec.binder_mass_c = (2 * mass_cathode * Pouch.jr_area) * Binder_share # 5% of total share
    getMass_elements_pouch_c_spec.conductive_material = (2 * mass_cathode * Pouch.jr_area) * Conductive_material_share # 5% of total share
    getMass_elements_pouch_c_spec.Separator_mass = 2 * getMass_separator.mass_sep * Pouch.jr_area
    
    return getMass_elements_pouch_c_spec.Ni_mass, getMass_elements_pouch_c_spec.Co_mass, \
           getMass_elements_pouch_c_spec.Mn_mass, getMass_elements_pouch_c_spec.Li_mass, \
           getMass_elements_pouch_c_spec.Al_mass

def getMass_elements_pouch_a(Electrodecomposition_anode_opt1, CurrentCollector_anode, Pouch): 
    #
    mass_anode = Electrodecomposition_anode_opt1.active_load / Electrodecomposition_anode_opt1.active_frac
    mass_cc_a = CurrentCollector_anode.mass_cc_a
    # Graphite
    getMass_elements_pouch_a.Gr_mass = (2 * mass_anode * Pouch.jr_area) * Electrodecomposition_anode_opt1.active_frac
    getMass_elements_pouch_a.cc_a_mass = mass_cc_a * Pouch.jr_area
    getMass_elements_pouch_a.binder_mass_a = (2 * mass_anode * Pouch.jr_area) * 0.03 # 3% of total share
    return getMass_elements_pouch_a.Gr_mass, getMass_elements_pouch_a.cc_a_mass, \
           getMass_elements_pouch_a.binder_mass_a

def getMass_elements_pouch_elyte(Electrolyte_liquid, Separator, Pouch):
    #
    getMass_electrolyte_liquid.mass_elyte = Electrolyte_liquid.density_elyte * Pouch.ecap_ratio * Pouch.capacity
    if Electrolyte_liquid.salt == "LiPF6":
        Li_mass_content_LiPF6 = 0.0457 # 4.57 weight% of Li content in LiPF6 / concentration in mol L-1
        getMass_elements_cyl_elyte.Li_mass = getMass_electrolyte_liquid.mass_elyte_salt * Li_mass_content_LiPF6
    # elif Electrolyte.salt == "LLZO":
    #     getMass_separator.mass_sep = Separator.thickness_sep * Separator.density_sep  # g / cm2
    #     Li_mass_content_LLZO = 0.0579  # 5.79 weight% of Li content in LLZO
    #     getMass_elements_cyl_elyte.Li_mass = (Electrolyte.concentration * getMass_electrolyte.mass_elyte + getMass_separator.mass_sep) \
    #                                          * Li_mass_content_LLZO
    else:
        getMass_elements_cyl_elyte.Li_mass = 1
    return getMass_elements_cyl_elyte.Li_mass

#######################
#######################
# Get costs for all components Pouch
def getCosts_pouch(Ni_raw , Co_raw, Mn_raw, Li_raw, Al_raw, Al_cc_raw, Gr_raw, Cu_raw, \
                 Binder_raw, Elyte_raw, Separator_raw, Pouch_packaging_raw, Conductive_raw, Electrolyte_liquid):
    # Cathode
    Ni_costs = getMass_elements_pouch_c.Ni_mass/1000 * Ni_raw # gram / 1000 * $ / kg = $ 
    Co_costs = getMass_elements_pouch_c.Co_mass/1000 * Co_raw
    Mn_costs = getMass_elements_pouch_c.Mn_mass/1000 * Mn_raw
    Li_costs = getMass_elements_pouch_c.Li_mass/1000 * Li_raw
    Al_costs = getMass_elements_pouch_c.Al_mass/1000 * Al_raw
    Al_cc_costs = getMass_elements_pouch_c.cc_c_mass/1000 * Al_cc_raw
    # Anode
    Gr_costs = getMass_elements_pouch_a.Gr_mass/1000 * Gr_raw
    Cu_costs = getMass_elements_pouch_a.cc_a_mass/1000 * Cu_raw
    # Binder for both
    Binder_costs = (getMass_elements_pouch_c.binder_mass_c + getMass_elements_pouch_a.binder_mass_a)/1000 * Binder_raw
    # Electrolyte
    if Electrolyte_liquid.salt == "LiPF6":
        Elyte_costs = getMass_electrolyte_liquid.mass_elyte / 1000 * Elyte_raw
    else:
        Elyte_costs = getMass_elements_elyte_solid.Li_mass_abs / 1000 * Li_raw
    # Separator
    Separator_costs = getMass_elements_pouch_c.Separator_mass/1000 * Separator_raw
    # Housing
    Housing_costs = getMass_Pouch_total.mass_packaging_pouch/1000 * Pouch_packaging_raw
    # Conductive Material
    Conductive_material_costs = getMass_elements_pouch_c.conductive_material/1000 * Conductive_raw
    return Ni_costs, Co_costs, Mn_costs, Li_costs, Al_costs, Al_cc_costs, Gr_costs, Cu_costs, \
        Binder_costs, Elyte_costs, Separator_costs, Housing_costs, Conductive_material_costs

############################
# Mass for prismatic cell     
def getMass_elements_prismatic_c(Electrodecomposition_cathode_opt1, CurrentCollector_cathode, Activematerial_cathode, Prismatic):
    #
    mass_cathode = Electrodecomposition_cathode_opt1.active_load / Electrodecomposition_cathode_opt1.active_frac
    mass_cc_c = CurrentCollector_cathode.mass_cc_c
    # NMC811
    if Activematerial_cathode.name == "NMC811": # energies-434768-supplementary.xlsx
        getMass_elements_prismatic_c.Ni_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.48 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_prismatic_c.Co_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.06 * 0.9
        getMass_elements_prismatic_c.Mn_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.06 * 0.9
        getMass_elements_prismatic_c.Li_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.07 * 0.9
        getMass_elements_prismatic_c.Al_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.00 * 0.9
        getMass_elements_prismatic_c.cc_c_mass = mass_cc_c * Prismatic.jr_area
        getMass_elements_prismatic_c.binder_mass_c = (2 * mass_cathode * Prismatic.jr_area) * 0.05 # 5% of total share
        getMass_elements_prismatic_c.conductive_material = (2 * mass_cathode * Prismatic.jr_area) * 0.05 # 5% of total share
        getMass_elements_prismatic_c.Separator_mass = 2 * getMass_separator.mass_sep * Prismatic.jr_area
    # NMC622
    elif Activematerial_cathode.name == "NMC622":
        getMass_elements_prismatic_c.Ni_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.36 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_prismatic_c.Co_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.12 * 0.9
        getMass_elements_prismatic_c.Mn_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.12 * 0.9
        getMass_elements_prismatic_c.Li_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.07 * 0.9
        getMass_elements_prismatic_c.Al_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.00 * 0.9
        getMass_elements_prismatic_c.cc_c_mass = mass_cc_c * Prismatic.jr_area
        getMass_elements_prismatic_c.binder_mass_c = (2 * mass_cathode * Prismatic.jr_area) * 0.05 # 5% of total share
        getMass_elements_prismatic_c.conductive_material = (2 * mass_cathode * Prismatic.jr_area) * 0.05 # 5% of total share
        getMass_elements_prismatic_c.Separator_mass = 2 * getMass_separator.mass_sep * Prismatic.jr_area
    # NMC333
    elif Activematerial_cathode.name == "NMC333":
        getMass_elements_prismatic_c.Ni_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.20 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_prismatic_c.Co_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.20 * 0.9
        getMass_elements_prismatic_c.Mn_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.20 * 0.9
        getMass_elements_prismatic_c.Li_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.07 * 0.9
        getMass_elements_prismatic_c.Al_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.00 * 0.9
        getMass_elements_prismatic_c.cc_c_mass = mass_cc_c * Prismatic.jr_area
        getMass_elements_prismatic_c.binder_mass_c = (2 * mass_cathode * Prismatic.jr_area) * 0.05 # 5% of total share
        getMass_elements_prismatic_c.conductive_material = (2 * mass_cathode * Prismatic.jr_area) * 0.05 # 5% of total share
        getMass_elements_prismatic_c.Separator_mass = 2 * getMass_separator.mass_sep * Prismatic.jr_area
    # NCA
    elif Activematerial_cathode.name == "NCA":
        getMass_elements_prismatic_c.Ni_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.47 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_prismatic_c.Co_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.09 * 0.9
        getMass_elements_prismatic_c.Mn_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.00 * 0.9
        getMass_elements_prismatic_c.Li_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.07 * 0.9
        getMass_elements_prismatic_c.Al_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.01 * 0.9
        getMass_elements_prismatic_c.cc_c_mass = mass_cc_c * Prismatic.jr_area
        getMass_elements_prismatic_c.binder_mass_c = (2 * mass_cathode * Prismatic.jr_area) * 0.05 # 5% of total share
        getMass_elements_prismatic_c.conductive_material = (2 * mass_cathode * Prismatic.jr_area) * 0.05 # 5% of total share
        getMass_elements_prismatic_c.Separator_mass = 2 * getMass_separator.mass_sep * Prismatic.jr_area
    # LFP
    elif Activematerial_cathode.name == "LFP":
        getMass_elements_prismatic_c.Ni_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.00 * 0.9 # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
        getMass_elements_prismatic_c.Co_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.00 * 0.9
        getMass_elements_prismatic_c.Mn_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.00 * 0.9
        getMass_elements_prismatic_c.Li_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.04 * 0.9
        getMass_elements_prismatic_c.Al_mass = (2 * mass_cathode * Prismatic.jr_area) * 0.00 * 0.9
        getMass_elements_prismatic_c.cc_c_mass = mass_cc_c * Prismatic.jr_area
        getMass_elements_prismatic_c.binder_mass_c = (2 * mass_cathode * Prismatic.jr_area) * 0.05 # 5% of total share
        getMass_elements_prismatic_c.conductive_material = (2 * mass_cathode * Prismatic.jr_area) * 0.05 # 5% of total share
        getMass_elements_prismatic_c.Separator_mass = 2 * getMass_separator.mass_sep * Prismatic.jr_area
    else:
        1
    return getMass_elements_prismatic_c.Ni_mass, getMass_elements_prismatic_c.Co_mass, \
           getMass_elements_prismatic_c.Mn_mass, getMass_elements_prismatic_c.Li_mass, \
           getMass_elements_prismatic_c.Al_mass

# individual cell chemistries
def getMass_elements_prismatic_c_spec(Electrodecomposition_cathode_opt1, CurrentCollector_cathode, Activematerial_cathode, Prismatic,\
                                Ni_content, Co_content, Mn_content, Li_content, Al_content, Active_material_share, Binder_share, Conductive_material_share):
    #
    mass_cathode = Electrodecomposition_cathode_opt1.active_load / Electrodecomposition_cathode_opt1.active_frac
    mass_cc_c = CurrentCollector_cathode.mass_cc_c
    #
    getMass_elements_prismatic_c_spec.Ni_mass = (2 * mass_cathode * Prismatic.jr_area) * Ni_content * Active_material_share # twice the mass of cathode, due to double side coating, and 90% of cathode is made from AM
    getMass_elements_prismatic_c_spec.Co_mass = (2 * mass_cathode * Prismatic.jr_area) * Co_content * Active_material_share
    getMass_elements_prismatic_c_spec.Mn_mass = (2 * mass_cathode * Prismatic.jr_area) * Mn_content * Active_material_share
    getMass_elements_prismatic_c_spec.Li_mass = (2 * mass_cathode * Prismatic.jr_area) * Li_content * Active_material_share
    getMass_elements_prismatic_c_spec.Al_mass = (2 * mass_cathode * Prismatic.jr_area) * Al_content * Active_material_share
    getMass_elements_prismatic_c_spec.cc_c_mass = mass_cc_c * Pouch.jr_area
    getMass_elements_prismatic_c_spec.binder_mass_c = (2 * mass_cathode * Pouch.jr_area) * Binder_share # 5% of total share
    getMass_elements_prismatic_c_spec.conductive_material = (2 * mass_cathode * Pouch.jr_area) * Conductive_material_share # 5% of total share
    getMass_elements_prismatic_c_spec.Separator_mass = 2 * getMass_separator.mass_sep * Pouch.jr_area
    
    return getMass_elements_prismatic_c_spec.Ni_mass, getMass_elements_prismatic_c_spec.Co_mass, \
           getMass_elements_prismatic_c_spec.Mn_mass, getMass_elements_prismatic_c_spec.Li_mass, \
           getMass_elements_prismatic_c_spec.Al_mass

def getMass_elements_pristmatic_elyte(Electrolyte_liquid, Separator, Prismatic):
    #
    getMass_electrolyte_liquid.mass_elyte = Electrolyte_liquid.density_elyte * Prismatic.ecap_ratio * Prismatic.capacity
    if Electrolyte_liquid.salt == "LiPF6":
        Li_mass_content_LiPF6 = 0.0457 # 4.57 weight% of Li content in LiPF6 / concentration in mol L-1
        getMass_elements_cyl_elyte.Li_mass = getMass_electrolyte_liquid.mass_elyte_salt * Li_mass_content_LiPF6
    # elif Electrolyte.salt == "LLZO":
    #     getMass_separator.mass_sep = Separator.thickness_sep * Separator.density_sep  # g / cm2
    #     Li_mass_content_LLZO = 0.0579  # 5.79 weight% of Li content in LLZO
    #     getMass_elements_cyl_elyte.Li_mass = (Electrolyte.concentration * getMass_electrolyte.mass_elyte + getMass_separator.mass_sep) \
    #                                          * Li_mass_content_LLZO
    else:
        getMass_elements_cyl_elyte.Li_mass = 1
    return getMass_elements_cyl_elyte.Li_mass

#######################
#######################
# Get costs for all components Pouch
def getCosts_prismatic(Ni_raw , Co_raw, Mn_raw, Li_raw, Al_raw, Al_cc_raw, Gr_raw, Cu_raw, \
                 Binder_raw, Elyte_raw, Separator_raw, Pouch_packaging_raw, Conductive_raw):
    # Cathode
    Ni_costs = getMass_elements_prismatic_c.Ni_mass/1000 * Ni_raw # gram / 1000 * $ / kg = $ 
    Co_costs = getMass_elements_prismatic_c.Co_mass/1000 * Co_raw
    Mn_costs = getMass_elements_prismatic_c.Mn_mass/1000 * Mn_raw
    Li_costs = getMass_elements_prismatic_c.Li_mass/1000 * Li_raw
    Al_costs = getMass_elements_prismatic_c.Al_mass/1000 * Al_raw
    Al_cc_costs = getMass_elements_prismatic_c.cc_c_mass/1000 * Al_cc_raw
    # Anode
    Gr_costs = getMass_elements_prismatic_a.Gr_mass/1000 * Gr_raw
    Cu_costs = getMass_elements_prismatic_a.cc_a_mass/1000 * Cu_raw
    # Binder for both
    Binder_costs = (getMass_elements_prismatic_c.binder_mass_c + getMass_elements_prismatic_a.binder_mass_a)/1000 * Binder_raw
    # Electrolyte
    if Electrolyte_liquid.salt == "LiPF6" or "LiTFSI":
        Elyte_costs = getMass_electrolyte_liquid.mass_elyte / 1000 * Elyte_raw
    else:
        Elyte_costs = getMass_elements_elyte_solid.Li_mass_abs / 1000 * Li_raw
    # Separator
    Separator_costs = getMass_elements_prismatic_c.Separator_mass/1000 * Separator_raw
    # Housing
    Housing_costs = getMass_prismatic_total.mass_packaging_prismatic/1000 * Pouch_packaging_raw
    # Conductive Material
    Conductive_material_costs = getMass_elements_prismatic_c.conductive_material/1000 * Conductive_raw
    return Ni_costs, Co_costs, Mn_costs, Li_costs, Al_costs, Al_cc_costs, Gr_costs, Cu_costs, \
        Binder_costs, Elyte_costs, Separator_costs, Housing_costs, Conductive_material_costs
