#!/usr/bin/env python3
"""Write the synthetic 2D gate decks (gate/decks/twod/*.inp and their side files).

They reach the 2D input paths upstream's own decks do not: every [2D_OPTIONS] keyword and value,
quads, initial velocity, node maps by index and by tag, every boundary-condition type, edge
conveyance, the surface quality and infiltration rows, an external mesh file (with its own
`;; UNITS:` header, its own sections, a missing and a malformed one, a name with spaces), the
SI-units header in the .inp, the [GW_*], [2D_AQUIFER*], [2D_SURFACE_OWNERSHIP] and surface
quality sections and their resolution, the open() warnings and aliases, the 2D branches of the
1D code open() reaches (a rain gage only the mesh reads, the virtual-junction coupling rule, heat
overrides on the mesh), and one malformed row per diagnostic. They are gate inputs: both engines
must agree on each (gate/ctxgate.sh open2d, and every 1D gate, where the 2D sections are
unknown). Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
Run from anywhere: python3 gate/decks/twod/gen.py
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

BASE_1D = """[TITLE]
synthetic 2D open gate: {title}

[OPTIONS]
FLOW_UNITS {units}
FLOW_ROUTING DYNWAVE
START_DATE 01/01/2024
START_TIME 00:00:00
END_DATE 01/01/2024
END_TIME 01:00:00
REPORT_STEP 00:05:00
WET_STEP 00:05:00
DRY_STEP 00:05:00
ROUTING_STEP 00:00:10

[JUNCTIONS]
J1 99.0 1.5 0 0 0
J2 98.5 1.5 0 0 0

[OUTFALLS]
O1 98.0 FREE NO

[CONDUITS]
C1 J1 J2 50 0.013 0 0 0 0
C2 J2 O1 50 0.013 0 0 0 0

[XSECTIONS]
C1 CIRCULAR 0.5 0 0 0 1
C2 CIRCULAR 0.5 0 0 0 1

[TIMESERIES]
TS1 01/01/2024 00:00 100.2
TS1 01/01/2024 01:00 100.4
FLOW1 01/01/2024 00:00 0.01
FLOW1 01/01/2024 01:00 0.02

[CURVES]
RC1 RATING 0 0
RC1        1 1

[POLLUTANTS]
TSS MG/L 0 0 0 0 0
"""

MESH = """[2D_VERTICES]
0 0 100.0 va
10 0 100.0 vb
10 10 100.5 vc
0 10 100.5 vd
20 0 100.0 5001
20 10 100.5 vf

[2D_TRIANGLES]
0 1 2 0.030 0.10 low
0 2 3 0.035 west
"""
QUADS = """[2D_QUADS]
1 4 5 2 0.040 0.0 east
"""


def deck(name, title, body, units='CMS', mesh=True, quads=True, extra_1d='', files=None):
    text = BASE_1D.format(title=title, units=units) + extra_1d + '\n'
    if mesh:
        text += MESH + '\n'
        if quads:
            text += QUADS + '\n'
    text += body.strip('\n') + '\n'
    with open(os.path.join(HERE, name + '.inp'), 'w') as f:
        f.write(text)
    for fname, content in (files or {}).items():
        with open(os.path.join(HERE, fname), 'w') as f:
            f.write(content)


OPTIONS_A = """[2D_OPTIONS]
MAX_TIMESTEP 5.0
DRY_DEPTH 0.002
COUPLING_CD 0.6
DISPERSION 0.01
COUPLING_SYNC 30
LIMITER_EPSILON 1e-5
FLUX_DH_EPS 0.005
CELL_CLOSURE VFR
FACE_RECONSTRUCTION VFR_FACE
VFR_MIN_WET_FRAC 0.05
RAINFALL_MODE NEAREST_NEIGHBOR
REPORT_2D NO
OUTPUT_FILE results 2d.h5
INTEGRATOR EXPLICIT
MOMENTUM_EQUATION FULL_SWE
CORIOLIS LATITUDE
CORIOLIS_LATITUDE 45.5
CORIOLIS_F 1e-4
SMAGORINSKY_CS 0.2
EDDY_VISCOSITY_BACKGROUND 0.01
MIXING_LENGTH_COEFFICIENT 0.4
TURBULENCE_MODEL SMAGORINSKY
FRONT_REBUILD YES
RECONSTRUCTION_ORDER 2
THETA 0.9
THETA_REFERENCE_STEP 2.5
CFL_NUMBER 0.5
H_MOVE 0.004
LTS_TIERS 6
FROUDE_MAX 2.0
ADVECTION YES
COUPLING_IN_FLOODING YES
COUPLING_LAW WEIR_ORIFICE
COUPLING_HEAD STENCIL_WET_MEAN
COUPLING_POND_MODE TRACKING
COUPLING_POND_RELAX 60
COUPLING_RIM_BAND 0.25
COUPLING_EQUALIZE YES
COUPLING_COEFFS RUBINATO2017
COUPLING_AREA AUTO
BACKEND CPU
OUTPUT_PRECISION F64
OUTPUT_COMPRESSION 9
REPORT_2D_VARIABLES DEPTH, velocity HEAT
REPORT_2D_SPECIES TSS __TEMPERATURE__
REPORT_2D_STEP 00:10:00
INFILTRATION YES
INFIL_STEP 0:01:30
INFIL_DEFAULT_METHOD MODIFIED_GREEN_AMPT
INFIL_DESTINATION AQUIFER_2D
EVAPORATION CLIMATE
TRANSPORT_POLLUTANTS NO
TRANSPORT_MSX OFF
TRANSPORT_AGE NO
TRANSPORT_TEMPERATURE OFF
GROUNDWATER YES
GW_ET boundary_et
HEAT_MODEL COUPLED
HEAT_SURFACE_EXCHANGE YES
HEAT_RADIATIVE_EXCHANGE NO
BED_HEAT_MODEL GROUND_COLUMN
MIN_TIMESTEP 0.1
MAX_CVODE_STEPS 500
INTEGRATOR CVODE
"""

OPTIONS_B = """[2D_OPTIONS]
CELL_CLOSURE flat
FACE_RECONSTRUCTION mean
RAINFALL_MODE SYSTEM
REPORT_2D 1
MOMENTUM_EQUATION DW
CORIOLIS CONSTANT_F
CORIOLIS_LATITUDE UNSET
CORIOLIS_F unset
TURBULENCE_MODEL MIXING_LENGTH
FRONT_REBUILD off
RECONSTRUCTION_ORDER 1
ADVECTION false
COUPLING_IN_FLOODING TRUE
COUPLING_LAW ORIFICE
COUPLING_HEAD VERTEX
COUPLING_POND_MODE CHAMBER
COUPLING_RIM_BAND DEFAULT
COUPLING_EQUALIZE NO
COUPLING_AREA DEFAULT
BACKEND sycl
OUTPUT_PRECISION SINGLE
REPORT_2D_VARIABLES MINIMAL
REPORT_2D_SPECIES ALL
REPORT_2D_STEP 300
INFILTRATION NO
INFIL_DEFAULT_METHOD CURVE_NUM
INFIL_DESTINATION SUBCATCH_AQUIFER
EVAPORATION FORCING
GROUNDWATER NO
GW_ET NONE
HEAT_MODEL TRANSPORT_ONLY
HEAT_SURFACE_EXCHANGE AUTO
HEAT_RADIATIVE_EXCHANGE YES
BED_HEAT_MODEL LOCAL_STORAGE
REL_TOLERANCE 1e-4
ABS_TOLERANCE 1e-6
MAX_KRYLOV_DIM 5
LINEAR_SOLVER GMRES
PRECONDITIONER NONE
JACOBIAN AUTO
ATOL_AREA_REF 1
COUPLING_INTERVAL 1
COUPLING_WINDOW 1
ACTIVE_SET YES
ACTIVE_SET_HALO 2
MOMENTUM YES
"""

OPTIONS_C = """[2D_OPTIONS]
RAINFALL_MODE NATURAL_NEIGHBOUR
REPORT_2D yes
MOMENTUM_EQUATION LI
CORIOLIS OFF
TURBULENCE_MODEL NONE
FRONT_REBUILD AUTO
COUPLING_HEAD CELL
COUPLING_POND_MODE FOOTPRINT
COUPLING_EQUALIZE DEFAULT
COUPLING_COEFFS rubinato2017_upper
BACKEND HIP
OUTPUT_PRECISION DOUBLE
REPORT_2D_VARIABLES ALL
INFILTRATION AUTO
INFIL_DEFAULT_METHOD MODIFIED_HORTON
INFIL_DESTINATION LOST
EVAPORATION NO
GROUNDWATER AUTO
GW_ET CAPILLARY_RISE
BED_HEAT_MODEL NONE
MOMENTUM_EQUATION FULL
ADVECTION ON
REPORT_2D_VARIABLES NONE
"""

OPTIONS_D = """[2D_OPTIONS]
RAINFALL_MODE NONE
MOMENTUM_EQUATION SWE
MOMENTUM_EQUATION DIFFUSIVE
MOMENTUM_EQUATION DIFFUSIVE_WAVE
BACKEND OMP
BACKEND CUDA
BACKEND AUTO
OUTPUT_PRECISION FLOAT32
OUTPUT_PRECISION FLOAT64
REPORT_2D_VARIABLES DEFAULT GROUNDWATER GW_DETAILED TURBULENCE BUILDUP CONTINUITY GRADIENTS COUPLING INFILTRATION RAINFALL SPECIES NODE_HEAD EDGE_FLUX ENVELOPES
INFIL_DEFAULT_METHOD NONE
INFIL_DEFAULT_METHOD HORTON
INFIL_DEFAULT_METHOD MOD_HORTON
INFIL_DEFAULT_METHOD GREEN_AMPT
INFIL_DEFAULT_METHOD MOD_GREEN_AMPT
INFIL_DEFAULT_METHOD CURVE_NUMBER
INFIL_DEFAULT_METHOD CONSTANT
GW_ET AUTO
GW_ET BOTH
EVAPORATION ON
EVAPORATION OFF
EVAPORATION YES
INFILTRATION ON
INFILTRATION OFF
FRONT_REBUILD ON
FRONT_REBUILD NO
ADVECTION OFF
ADVECTION NO
ADVECTION TRUE
COUPLING_IN_FLOODING NO
COUPLING_IN_FLOODING FALSE
REPORT_2D 0
REPORT_2D NO
TRANSPORT_POLLUTANTS ON
TRANSPORT_MSX YES
TRANSPORT_AGE ON
TRANSPORT_TEMPERATURE YES
GROUNDWATER ON
GROUNDWATER OFF
MAX_TIMESTEP  +12.5
DRY_DEPTH 0x1p-8
COUPLING_CD inf
LIMITER_EPSILON nan
FLUX_DH_EPS   1e-400
H_MOVE 0
REPORT_2D_STEP 1:30
INFIL_STEP 45
LTS_TIERS +3
RECONSTRUCTION_ORDER 0002
OUTPUT_COMPRESSION 0
"""

# Each bad [2D_OPTIONS] value: its own deck, so each diagnostic is compared.
BAD_OPTIONS = [
    'MAX_TIMESTEP x', 'MAX_TIMESTEP 5s', 'DRY_DEPTH', 'COUPLING_CD abc', 'DISPERSION -1',
    'DISPERSION inf', 'DISPERSION x', 'COUPLING_SYNC -2', 'COUPLING_SYNC q', 'LIMITER_EPSILON e',
    'FLUX_DH_EPS 1e', 'CELL_CLOSURE ROUND', 'FACE_RECONSTRUCTION MAX', 'VFR_MIN_WET_FRAC 0',
    'VFR_MIN_WET_FRAC 0.6', 'VFR_MIN_WET_FRAC z', 'RAINFALL_MODE IDW', 'REPORT_2D MAYBE',
    'INTEGRATOR ARKODE', 'MOMENTUM_EQUATION KINEMATIC', 'CORIOLIS ON',
    'CORIOLIS_LATITUDE 91', 'CORIOLIS_LATITUDE inf', 'CORIOLIS_LATITUDE x', 'CORIOLIS_F nan',
    'SMAGORINSKY_CS -0.1', 'SMAGORINSKY_CS UNSET', 'EDDY_VISCOSITY_BACKGROUND -1',
    'MIXING_LENGTH_COEFFICIENT -1', 'TURBULENCE_MODEL KEPS', 'FRONT_REBUILD MAYBE',
    'RECONSTRUCTION_ORDER 3', 'RECONSTRUCTION_ORDER 1.0', 'RECONSTRUCTION_ORDER x',
    'THETA 0', 'THETA 1.5', 'THETA nan', 'THETA_REFERENCE_STEP -1', 'THETA_REFERENCE_STEP inf',
    'CFL_NUMBER 0', 'CFL_NUMBER 2', 'H_MOVE -0.1', 'LTS_TIERS 0', 'LTS_TIERS 9', 'LTS_TIERS 4294967300',
    'FROUDE_MAX 0', 'ADVECTION MAYBE', 'COUPLING_IN_FLOODING ON', 'COUPLING_LAW WEIR',
    'COUPLING_HEAD NODE', 'COUPLING_POND_MODE SPONGE', 'COUPLING_POND_RELAX -5',
    'COUPLING_RIM_BAND -1', 'COUPLING_EQUALIZE ON', 'COUPLING_COEFFS FANCY', 'COUPLING_AREA BIG',
    'BACKEND GPU', 'OUTPUT_PRECISION HALF', 'OUTPUT_COMPRESSION 10', 'OUTPUT_COMPRESSION -1',
    'REPORT_2D_VARIABLES DEPTH SPEED', 'REPORT_2D_STEP -5', 'REPORT_2D_STEP inf',
    'INFILTRATION MAYBE', 'INFIL_STEP -3', 'INFIL_DEFAULT_METHOD PHILIP',
    'INFIL_DESTINATION SEWER', 'HEAT_MODEL FULL', 'HEAT_SURFACE_EXCHANGE MAYBE',
    'HEAT_RADIATIVE_EXCHANGE 1', 'BED_HEAT_MODEL SAND', 'EVAPORATION MAYBE',
    'TRANSPORT_POLLUTANTS MAYBE', 'TRANSPORT_AGE 1', 'GROUNDWATER MAYBE', 'GW_ET SOME',
    'NOT_A_KEY 1', 'MAX_TIMESTEP',
]

MESH_ROWS = """[2D_INITIAL_VELOCITY]
0 0.1 -0.2
2 0x10 1e-3

[2D_VERTEX_NODE_MAP]
2 J1 0.6 2.5
5001 J2
vd J2 x 3.0
vb O1 0.5 y

[2D_TRIANGLE_NODE_MAP]
low J1 0.55 4.0
1 J2
low J2 z w
east O1 0.7

[2D_BOUNDARY_CONDITIONS]
0 2 NORMAL_FLOW 0.005 * *
1 1 WALL 123 x wallgrp
1 2 WALL_DRAG 0.02
2 0 TS_STAGE TS1 * tide
2 1 SPECIFIED_STAGE 100.5
2 2 SPECIFIED_FLOW 0.05 * up
2 3 TS_FLOW FLOW1
0 0 RATING_CURVE RC1 * rc
1 0 NORMAL_FLOW *
0 1 SPECIFIED_STAGE

[2D_EDGE_CONVEYANCE]
0 2 0.5
2 1 1
4 5 0

[2D_INITIAL_QUALITY]
* TSS 1.5
CELL 2 TSS 3
TAG west __TEMPERATURE__ -2.5
TAG east __temperature__ -1
CELL 3 __WATER_AGE__ 4

[2D_INITIAL_BED_QUALITY]
* TSS 0.5
CELL 1 TSS 2
TAG low TSS 0

[2D_BOUNDARY_QUALITY]
2 0 TSS 10
2 3 __TEMPERATURE__ -4
0 2 TSS 0
"""

INFIL_ROWS = """[2D_OPTIONS]
INFIL_STEP 0:02:00
INFIL_DEFAULT_METHOD HORTON

[2D_INFILTRATION_OPTIONS]
INFIL_STEP 0:05

[2D_INFILTRATION_DEFAULTS]
* HORTON 3.0 0.5 4 7 -
low GREEN_AMPT 3.5 0.5 0.25 SUBCATCH_AQUIFER
west MODIFIED_HORTON 1 2
east NONE
* CONSTANT 0.5 AQUIFER_2D

[2D_INFILTRATION]
1 CURVE_NUMBER 80 - 7 LOST
3 MODIFIED_GREEN_AMPT 4 0.4 0.3
2 CONSTANT 1.25
"""

BAD_ROWS = {
    '2D_VERTICES': ['1 2', 'x 0 0', '0 y 0', '0 0 z', '0 0 1e'],
    '2D_TRIANGLES': ['0 1', 'a 1 2 0.03', '0 b 2 0.03', '0 1 c 0.03', '0 1 2 n', '0 1 2 0.03 -1', '0 1.5 2 0.03'],
    '2D_QUADS': ['0 1 2 3', 'a 1 2 3 0.03', '0 1 2 d 0.03', '0 1 2 3 n', '0 1 2 3 0.03 -0.5'],
    '2D_INITIAL_VELOCITY': ['0 1', '9 0 0', '-1 0 0', 'x 0 0', '0 u 0', '0 0 v'],
    '2D_VERTEX_NODE_MAP': ['0', 'nosuchtag J1', '99 J1', '-1 J1'],
    '2D_TRIANGLE_NODE_MAP': ['0', 'nosuch J1', '99 J1'],
    '2D_BOUNDARY_CONDITIONS': ['0 1', '-1 0 WALL', 'x 0 WALL', '0 4 WALL', '0 -1 WALL', '0 1 OPEN',
                               '0 1 WALL_DRAG', '0 1 WALL_DRAG -1', '0 1 WALL_DRAG inf', '0 1 NORMAL_FLOW x',
                               '0 1 SPECIFIED_STAGE 1e'],
    '2D_EDGE_CONVEYANCE': ['0 1', '-1 2 0.5', '0 x 0.5', '1 1 0.5', '0 1 c', '0 1 1.5', '0 1 -0.1'],
    '2D_INITIAL_QUALITY': ['* TSS', 'CELL 0 TSS 1', 'CELL 1 TSS', 'TAG low TSS', 'ALL TSS 1', '* TSS -1',
                           '* TSS inf', 'CELL x TSS 1', '* TSS c'],
    '2D_INITIAL_BED_QUALITY': ['* TSS', 'CELL 0 TSS 1', 'CELL 1 TSS', 'TAG low TSS', 'NONE TSS 1',
                               '* __TEMPERATURE__ -1', '* TSS nan'],
    '2D_BOUNDARY_QUALITY': ['0 1 TSS', '-1 0 TSS 1', '0 4 TSS 1', '0 1 TSS -2', '0 1 TSS x'],
    '2D_INFILTRATION_OPTIONS': ['INFIL_STEP', 'STEP 1', 'INFIL_STEP -1'],
    '2D_INFILTRATION_DEFAULTS': ['*', '* PHILIP 1', '* HORTON 1 2 3 4 5 6', '* HORTON 1 x LOST',
                                 '* HORTON 1 2 SEWER', 'low MOD_HORTON 1'],
    '2D_INFILTRATION': ['1', '0 HORTON 1', 'x HORTON 1', '1 HORTON 1 2 3 4 5 6 7', '1 GREEN_AMPT 1 q'],
}

GW_ROWS = """[GW_TRANSPORT_OPTIONS]
UNSAT_TRANSPORT layered
THERMAL_MODEL POROUS_EQUILIBRIUM
THERMAL_LAYERS 32
LATERAL_CONDUCTION ON
SOIL_ET_HEAT SENSIBLE_ONLY
THERMAL_MAX_STEP 30
WATER_THERMAL_CONDUCTIVITY 0.6
AIR_THERMAL_CONDUCTIVITY 0.03
SURFACE_THERMAL_RESISTANCE 0.2
TRANSPORT_POLLUTANTS NO
TRANSPORT_MSX OFF
TRANSPORT_AGE YES
TRANSPORT_TEMPERATURE ON
DISPERSION NO
CONDUCTION YES
SURFACE_THERMAL_BC FIXED 12.5
DEEP_THERMAL_BC FIXED_TEMP 10 DEPTH 5
THERMAL_MIXING GEOMETRIC
C_DIFF 0.5

[GW_TRANSPORT_PARAMS]
* 2650 880 2 0 0.5 0.05 1e-9 1e-9 0.065
TAG low - - 1.5
CELL 2 2000 - - - - - - - -12

[GW_SORPTION]
* TSS 0.5
TAG west TSS 0 0.1
CELL 3 TSS 1.5 -

[GW_INITIAL_QUALITY]
* SAT TSS 5
TAG low UNSAT TSS 2
CELL 2 LAYER 3 __TEMPERATURE__ -1
CELL 1 SAT __WATER_AGE__ -2
FILE tw_gw_iq.csv

[GW_BOUNDARY_QUALITY]
1 0 TSS CONC 1.5
2 3 TSS TS TS1
3 1 __TEMPERATURE__ HEATFLUX 1e3
1 2 TSS MASSFLUX someseries

[GW_SOURCES]
well1 CELL 2 FLOW 0.01 SCALE 2 TSS CONC 5 __TEMPERATURE__ MASS FLOW1
well2 XY 15 5 FLOW FLOW1
well3 * FLOW -0.001
well4 GLOBAL FLOW 0
well5 TAG low FLOW 1e-3 TSS MASS 0.5
"""
GW_IQ_CSV = """element,zone,species,value
CELL 1,SAT,TSS,3.5
TAG low;UNSAT;TSS;1
* LAYER 2 TSS 4
CELL 9 SAT TSS 1
* SAT NOPE 1
bad line here
CELL 1 SAT TSS -3
"""

GW_RESOLVE_ERRS = """[GW_TRANSPORT_PARAMS]
CELL 9 2650
TAG nosuchtag 2650

[GW_SORPTION]
* NOPE 0.5

[GW_INITIAL_QUALITY]
* SAT TSS -1
* SAT __TEMPERATURE__ -1
CELL 99 SAT TSS 1
FILE missing_sidecar.csv

[GW_BOUNDARY_QUALITY]
1 3 TSS CONC 1
3 3 TSS CONC 1
1 0 NOPE TS NOSERIES

[GW_SOURCES]
dup CELL 1 FLOW 0.01
dup CELL 2 FLOW NOFLOW
xy XY 1e3 -5 FLOW 0
"""

GW_BAD = {
    'GW_TRANSPORT_OPTIONS': ['UNSAT_TRANSPORT', 'UNSAT_TRANSPORT BULK x', 'UNSAT_TRANSPORT MIXED',
                             'THERMAL_MODEL HOT', 'THERMAL_LAYERS 0', 'THERMAL_LAYERS 2.5', 'THERMAL_LAYERS 129',
                             'LATERAL_CONDUCTION MAYBE', 'SOIL_ET_HEAT NONE', 'THERMAL_MAX_STEP 0',
                             'AIR_THERMAL_CONDUCTIVITY nan', 'TRANSPORT_AGE MAYBE', 'CONDUCTION 1',
                             'SURFACE_THERMAL_BC FIXED', 'SURFACE_THERMAL_BC FIXED hot', 'SURFACE_THERMAL_BC SUN',
                             'DEEP_THERMAL_BC MANTLE', 'DEEP_THERMAL_BC FIXED_TEMP', 'DEEP_THERMAL_BC FIXED_TEMP warm',
                             'DEEP_THERMAL_BC GEOTHERMAL_FLUX 0.06 DEPTH -1', 'THERMAL_MIXING HARMONIC',
                             'C_DIFF 0', 'C_DIFF 1.5', 'NOT_A_KEY 1'],
    'GW_TRANSPORT_PARAMS': ['ALL 1', 'TAG', 'CELL', 'CELL 0 1', 'CELL 1.5 1', '* 1 2 3 4 5 6 7 8 9 10',
                            '* x', '* -1'],
    'GW_SORPTION': ['* TSS', '* TSS -1', '* TSS 1 -1', '* TSS nan'],
    'GW_INITIAL_QUALITY': ['FILE', '* SAT TSS', '* WET TSS 1', '* LAYER 1 TSS', '* LAYER 0 TSS 1',
                           '* SAT TSS x'],
    'GW_BOUNDARY_QUALITY': ['1 0 TSS CONC', '0 0 TSS CONC 1', '1 -1 TSS CONC 1', '1 0 TSS FLUX 1'],
    'GW_SOURCES': ['w CELL 1', 'w CELL 0 FLOW 1', 'w XY 1 FLOW', 'w XY a b FLOW 1', 'w HERE FLOW 1',
                   'w * RATE 1 x', 'w * FLOW', 'w * FLOW 1 SCALE -1', 'w * FLOW 1 TSS CONC',
                   'w * FLOW 1 TSS LOAD 1'],
}

AQ_ROWS = """[2D_AQUIFER_OPTIONS]
GW_CLOCK INDEPENDENT
GW_BASE_STEP 2
GW_MAX_STEP 600
GW_BINDS_SURFACE_DT YES
GW_BUDGET_FRACTION 1
GW_HEAD_TOLERANCE 0.1
SOIL_CHAR brooks-corey
CLOSURE SIGMA
M_LAYERS 16
CAPILLARY_DIFF TRUE
FORCE_CLOSED_FORM 0
DUNNE NO
C_GW 0.4
C_COL 1
MODE PER_SUBCATCH
NODE_ENROLMENT ROWS
LINK_SEEPAGE ONE_WAY
GW_ET capillary_rise
WILTING_SUCTION 120

[2D_AQUIFER]
* 1e-5 5 0.45 0.1 2 PSI_B 0.3 LAMBDA 0.5 N 1.8 L 0.4 C_LOSS 0.01 HG0 2 SOIL_CHAR VG CLOSURE A M_LAYERS 8
TAG low 2e-5 4 0.4 0.05 1.5
CELL 3 1e-6 3 0.35 0 1 SOIL_CHAR Gardner CLOSURE ENSLAVED

[2D_AQUIFER_NODE]
J1 AUTO KC 1e-6 DC 0.5 AREA 2
J2 2 EXCHANGE NO
O1 1

[2D_AQUIFER_LINKS]
C1 KC 1e-7 DC 0.3
C2 EXCHANGE FALSE

[2D_AQUIFER_BOUNDARIES]
CELL 1 EDGE 0 HEAD CONST 99.5
CELL 2 EDGE 3 FLOW SERIES FLOW1

[2D_SURFACE_OWNERSHIP]
S1 SUBCATCH UNIFORM
s2 subcatch uniform
"""

AQ_BAD = {
    '2D_AQUIFER_OPTIONS': ['GW_CLOCK', 'GW_CLOCK SOMETIMES', 'GW_BASE_STEP 0', 'GW_MAX_STEP x',
                           'GW_BINDS_SURFACE_DT 2', 'GW_BUDGET_FRACTION 0', 'GW_HEAD_TOLERANCE 1',
                           'SOIL_CHAR CLAY', 'CLOSURE OPEN', 'M_LAYERS 1', 'M_LAYERS 4.5',
                           'CAPILLARY_DIFF Y', 'FORCE_CLOSED_FORM N', 'DUNNE MAYBE', 'C_GW 0',
                           'C_COL 2', 'MODE GRID', 'NODE_ENROLMENT ALL', 'LINK_SEEPAGE BOTH',
                           'GW_ET SOME', 'WILTING_SUCTION 0', 'WILTING_SUCTION inf', 'NOT_A_KEY 1'],
    '2D_AQUIFER': ['ALL 1 1 1 1 1', 'TAG', 'CELL 0 1 1 .4 .1 1', '* 1 1 .4 .1', '* x 1 .4 .1 1',
                   '* 1 1 .4 .1 1 PSI_B', '* 1 1 .4 .1 1 FOO 1', '* 1 1 .4 .1 1 SOIL_CHAR X',
                   '* 1 1 .4 .1 1 CLOSURE X', '* 1 1 .4 .1 1 M_LAYERS 200', '* 0 1 .4 .1 1',
                   '* 1 0 .4 .1 1', '* 1 1 1.4 .1 1', '* 1 1 .4 .5 1', '* 1 1 .4 .1 0', '* 1 1 .4 .1 1 N 1',
                   '* 1 1 .4 .1 1 LAMBDA 0', '* 1 1 .4 .1 1 PSI_B -1', '* 1 1 .4 .1 1 C_LOSS -1',
                   '* 1 1 .4 .1 1 N x', '* 1 1 .4 .1 1 HG0 y'],
    '2D_AQUIFER_NODE': ['J1', 'J1 0', 'J1 x', 'J1 1 KC', 'J1 1 KC x', 'J1 1 DC x', 'J1 1 AREA x',
                        'J1 1 EXCHANGE MAYBE', 'J1 1 FOO 1', 'J1 1 KC -1', 'J1 1 KC 1'],
    '2D_AQUIFER_LINKS': ['C1 KC', 'C1 KC x', 'C1 DC x', 'C1 EXCHANGE 2', 'C1 FOO 1', 'C1 DC -1'],
    '2D_AQUIFER_BOUNDARIES': ['CELL 1 EDGE 0 HEAD CONST', 'NODE 1 EDGE 0 HEAD CONST 1',
                              'CELL 0 EDGE 0 HEAD CONST 1', 'CELL 1 EDGE 4 HEAD CONST 1',
                              'CELL 1 EDGE 0 LEVEL CONST 1', 'CELL 1 EDGE 0 HEAD CONST inf',
                              'CELL 1 EDGE 0 HEAD SERIES ""', 'CELL 1 EDGE 0 HEAD TABLE x',
                              'CELL 1 EDGE 0 HEAD CONST 1\nCELL 1 EDGE 0 FLOW CONST 2'],
    '2D_SURFACE_OWNERSHIP': ['S1 SUBCATCH', 'S1 MESH UNIFORM', 'S1 SUBCATCH WEIGHTED',
                             'S1 SUBCATCH UNIFORM\ns1 SUBCATCH UNIFORM'],
}

SQ_ROWS = """[2D_COVERAGES]
* Res 60 Com 40
TAG low Res 100
CELL 2 Com 30
* Res 50 Com 50

[2D_LOADINGS]
* TSS 1.5
CELL 1 TSS 0.25
TAG west TSS 0

[2D_CURB_LENGTH]
* 120
CELL 3 40.5
"""
SQ_BAD = {
    '2D_COVERAGES': ['*', '* Res', '* Res x', '* Res -5', '* Res 60 Com 50', 'TAG', 'CELL', 'CELL 0 Res 1',
                     'CELL 1.5 Res 1', 'ALL Res 1', '* Res 1e400'],
    '2D_LOADINGS': ['* TSS', '* TSS 1 2', '* TSS -1', '* TSS 1e-320', 'CELL x TSS 1'],
    '2D_CURB_LENGTH': ['*', '* 1 2', '* -1', '* x'],
}

MESHFILE_2DM = """;; external mesh: its own units header
;; UNITS: SI (m)
[2D_OPTIONS]
DRY_DEPTH 0.005
MIN_TIMESTEP 0.5
[2D_VERTICES]
0 0 50 sa
5 0 50 sb
5 5 51 sc
0 5 51 sd
[2D_TRIANGLES]
0 1 2 0.02 0.05 side
0 2 3 0.02
[2D_QUADS]
1 0 3 2 0.03
[2D_INITIAL_VELOCITY]
1 0.2 0.1
[2D_VERTEX_NODE_MAP]
sc J1 0.7 1.5
[2D_TRIANGLE_NODE_MAP]
side J2 0.6
[2D_BOUNDARY_CONDITIONS]
0 1 SPECIFIED_STAGE 50.2
[2D_EDGE_CONVEYANCE]
0 2 0.25
[2D_INITIAL_QUALITY]
CELL 1 TSS 7
[2D_INITIAL_BED_QUALITY]
* TSS 1
[2D_BOUNDARY_QUALITY]
0 1 TSS 2
[2D_INFILTRATION_DEFAULTS]
* CONSTANT 0.75
[2D_INFILTRATION]
2 HORTON 1 0.5 2
[2D_INFILTRATION_OPTIONS]
INFIL_STEP 0:00:30
[JUNCTIONS]
IGNORED 1 1 0 0 0
[2D_MESH_FILE]
FILE recursion.2dm
"""
MESHFILE_INLINE = """[2D_OPTIONS]
RAINFALL_MODE SYSTEM
INFIL_STEP 0:10:00

[2D_MESH_FILE]
FILE tw_meshfile.2dm
FILE second.2dm

[2D_INFILTRATION_DEFAULTS]
* HORTON 3 0.5 4 7 0
west GREEN_AMPT 3 0.4 0.2

[2D_INITIAL_QUALITY]
* TSS 0.5
"""
MESHFILE_BAD_2DM = """[2D_VERTICES]
0 0 1
1 0 1
[2D_TRIANGLES]
0 1 2 0.03
0 1 bad 0.03
[2D_INFILTRATION_DEFAULTS]
* HORTON 1
"""


def main():
    # every [2D_OPTIONS] keyword, valid values (four decks: the enum variants and the warnings)
    deck('tw_options_a', 'every [2D_OPTIONS] key, first values; retired keys warn', OPTIONS_A)
    deck('tw_options_b', '[2D_OPTIONS] second values, the rest of the retired keys', OPTIONS_B)
    deck('tw_options_c', '[2D_OPTIONS] third values', OPTIONS_C)
    deck('tw_options_d', '[2D_OPTIONS] the remaining spellings and strtod/strtol edges', OPTIONS_D)
    for i, row in enumerate(BAD_OPTIONS):
        deck(f'tw_badopt_{i:02d}', f'bad [2D_OPTIONS] row: {row}', '[2D_OPTIONS]\nDRY_DEPTH 0.002\n' + row + '\nMAX_TIMESTEP 7\n')
    # the mesh and its row sections, every valid form
    deck('tw_mesh_rows', 'quads, initial velocity, node maps, boundary rows, conveyance, surface quality rows', MESH_ROWS)
    deck('tw_mesh_si_header', 'the ;; UNITS: SI (m) header mid-file in the .inp', MESH_ROWS + '\n;; UNITS: metres\n', units='CFS')
    deck('tw_mesh_si_variants', 'UNITS markers that do and do not count', ';;units:  FT\n;; UNITS:SI (M)\n;UNITS: m\n', units='LPS')
    deck('tw_mesh_ft_header', 'a non-SI UNITS marker leaves the flag', ';; UNITS: ft\n;; Units: feet\n', units='CFS')
    deck('tw_tri_after_quad', 'a triangle after a quad', '[2D_TRIANGLES]\n1 2 3 0.03\n')
    deck('tw_no_mesh_rows', 'node map and boundary rows without a mesh', '[2D_VERTEX_NODE_MAP]\n0 J1\n', mesh=False)
    deck('tw_no_mesh_bc', 'boundary, conveyance and quality rows without a mesh', '[2D_BOUNDARY_CONDITIONS]\n0 0 WALL\n[2D_EDGE_CONVEYANCE]\n0 1 0.5\n[2D_INITIAL_QUALITY]\n* TSS 1\n', mesh=False)
    for sec, rows in BAD_ROWS.items():
        for i, row in enumerate(rows):
            deck(f'tw_bad_{sec[3:].lower()}_{i:02d}', f'bad [{sec}] row: {row}', f'[{sec}]\n{row}\n')
    # infiltration rows and the INFIL_STEP alias
    deck('tw_infil', 'infiltration options, defaults and per-cell rows; the INFIL_STEP alias wins', INFIL_ROWS)
    deck('tw_infil_noalias', 'the infiltration INFIL_STEP without the alias', INFIL_ROWS.replace('INFIL_STEP 0:02:00\n', ''))
    # the external mesh file
    deck('tw_meshfile', 'an external mesh file with its own UNITS header and sections', MESHFILE_INLINE,
         mesh=False, files={'tw_meshfile.2dm': MESHFILE_2DM})
    deck('tw_meshfile_inline_mesh', 'an external mesh file after an inline mesh', MESHFILE_INLINE.replace('[2D_INFILTRATION_DEFAULTS]\n* HORTON 3 0.5 4 7 0\nwest GREEN_AMPT 3 0.4 0.2\n', ''))
    deck('tw_meshfile_spaces', 'a mesh file name with spaces', '[2D_MESH_FILE]\nfile   tw  mesh   file.2dm\n', mesh=False,
         files={'tw  mesh   file.2dm': MESHFILE_2DM.replace('[2D_MESH_FILE]\nFILE recursion.2dm\n', '')})
    deck('tw_meshfile_quoted', 'a quoted mesh file name', '[2D_MESH_FILE]\nNOTFILE x\nFILE "tw  mesh   file.2dm"\n', mesh=False)
    deck('tw_meshfile_missing', 'a missing mesh file (strict fails, lenient goes on 1D-only)', '[2D_MESH_FILE]\nFILE nowhere/none.2dm\n[2D_INITIAL_QUALITY]\n* TSS 2\n[2D_BOUNDARY_CONDITIONS]\n0 0 WALL\n')
    deck('tw_meshfile_bad', 'a mesh file with a malformed row', '[2D_MESH_FILE]\nFILE tw_meshfile_bad.2dm\n', mesh=False,
         files={'tw_meshfile_bad.2dm': MESHFILE_BAD_2DM})
    deck('tw_meshfile_dir', 'a directory as the mesh file', '[2D_MESH_FILE]\nFILE .\n', mesh=False)
    deck('tw_meshfile_updir', 'a mesh file path through ..', '[2D_MESH_FILE]\nFILE ../twod/./tw_meshfile_bad.2dm\n', mesh=False)
    # groundwater transport rows, their resolution and the inert warning
    deck('tw_gw', 'the [GW_*] rows, a FILE of initial quality, XY sources; inert without [2D_AQUIFER]', GW_ROWS,
         files={'tw_gw_iq.csv': GW_IQ_CSV})
    deck('tw_gw_resolve_errs', 'every resolveGwTransport diagnostic', GW_RESOLVE_ERRS)
    deck('tw_gw_xy_nomesh', 'an XY source on a deck without cells', '[GW_SOURCES]\nw XY 1 2 FLOW 0\n', mesh=False)
    deck('tw_gw_opts_only', 'an options-only [GW_*] block counts as authored', '[GW_TRANSPORT_OPTIONS]\nSURFACE_THERMAL_BC ADIABATIC\nDEEP_THERMAL_BC ADIABATIC\n')
    deck('tw_gw_bcs', 'the thermal boundary spellings', '[GW_TRANSPORT_OPTIONS]\nSURFACE_THERMAL_BC SURFACE_WATER\nSURFACE_THERMAL_BC HEAT_FLUX 5\nSURFACE_THERMAL_BC TIMESERIES TS1\nSURFACE_THERMAL_BC FLUX_TIMESERIES TS1\nSURFACE_THERMAL_BC AIR\nDEEP_THERMAL_BC TIMESERIES TS1 DEPTH 2\nDEEP_THERMAL_BC FLUX_TIMESERIES TS1\nDEEP_THERMAL_BC GEOTHERMAL_FLUX 0.07 DEPTH\nTHERMAL_MIXING ARITHMETIC\nTHERMAL_MODEL WATER_ONLY\nUNSAT_TRANSPORT BULK\nLATERAL_CONDUCTION NO\nSOIL_ET_HEAT LATENT_SURFACE\n')
    for sec, rows in GW_BAD.items():
        for i, row in enumerate(rows):
            deck(f'tw_badgw_{sec[3:].lower()}_{i:02d}', f'bad [{sec}] row: {row}', f'[{sec}]\n{row}\n')
    # aquifer, ownership; the GW_ET alias; GROUNDWATER YES without rows
    deck('tw_aquifer', 'the [2D_AQUIFER*] and [2D_SURFACE_OWNERSHIP] rows', AQ_ROWS + '\n[2D_OPTIONS]\nGW_ET BOTH\n' + GW_ROWS.split('[GW_SORPTION]')[0])
    deck('tw_gw_et_alias', 'the [2D_OPTIONS] GW_ET alias alone authors the aquifer options', '[2D_OPTIONS]\nGW_ET NONE\n')
    deck('tw_groundwater_yes', 'GROUNDWATER YES with no aquifer rows warns', '[2D_OPTIONS]\nGROUNDWATER YES\n')
    for sec, rows in AQ_BAD.items():
        for i, row in enumerate(rows):
            deck(f'tw_badaq_{sec[3:].lower()}_{i:02d}', f'bad [{sec}] row', f'[{sec}]\n{row}\n')
    # surface quality rows
    deck('tw_sq', 'the [2D_COVERAGES] / [2D_LOADINGS] / [2D_CURB_LENGTH] rows', SQ_ROWS)
    for sec, rows in SQ_BAD.items():
        for i, row in enumerate(rows):
            deck(f'tw_badsq_{sec[3:].lower()}_{i:02d}', f'bad [{sec}] row: {row}', f'[{sec}]\n{row}\n')
    # 2D branches of the 1D code open() reaches
    rain = '[RAINGAGES]\nRG1 INTENSITY 0:05 1.0 TIMESERIES TS1\nRG2 INTENSITY 0:05 1.0 TIMESERIES FLOW1\n'
    deck('tw_rain_mesh', 'gages no subcatchment names are used under a rain-on-grid mesh', '[2D_OPTIONS]\nRAINFALL_MODE SYSTEM\n', extra_1d=rain)
    deck('tw_rain_none', 'RAINFALL_MODE NONE leaves unused gages unused', '[2D_OPTIONS]\nRAINFALL_MODE NONE\n', extra_1d=rain)
    deck('tw_rain_nomesh', 'no mesh: unused gages stay unused', '[2D_OPTIONS]\nRAINFALL_MODE SYSTEM\n', extra_1d=rain, mesh=False)
    vj = '[VIRTUAL_JUNCTIONS]\nV1 98.8\nV2 98.6\n\n[CONDUITS]\nK1 J1 V1 30 0.013 0 0 0 0\nK2 V1 J2 30 0.013 0 0 0 0\nK3 J1 V2 30 0.013 0 0 0 0\nK4 V2 J2 30 0.013 0 0 0 0\n\n[XSECTIONS]\nK1 CIRCULAR 0.5 0 0 0 1\nK2 CIRCULAR 0.5 0 0 0 1\nK3 CIRCULAR 0.5 0 0 0 1\nK4 CIRCULAR 0.5 0 0 0 1\n'
    deck('tw_vj_coupled', 'virtual junctions a node map couples to the mesh (rule 5)', '[2D_VERTEX_NODE_MAP]\n2 V1\n[2D_TRIANGLE_NODE_MAP]\nwest V2\n', extra_1d=vj)
    deck('tw_vj_uncoupled', 'virtual junctions with a mesh but no coupling', '[2D_VERTEX_NODE_MAP]\n2 J1\n', extra_1d=vj)
    heat = '[PROCESS_COMPONENTS]\norg.hydrocouple.openswmm.heat config="tw_heat.heat"\n\n'
    heat_cfg = '[RADIATIVE_FLUXES]\nALBEDO CELL2D 2 0.3\nSKY_VIEW MESH_TAG low 0.6\nSHADE_FACTOR CELL2D 0 0.5\nALBEDO CELL2D 1x 0.5\nALBEDO CELL2D 9 0.5\nALBEDO MESH_TAG nosuch 0.5\n[SEDIMENT_EXCHANGE]\nBED_THICKNESS MESH_TAG east 0.4\nBED_THICKNESS CELL2D 1 0.2\n'
    deck('tw_heat_mesh', 'heat overrides on mesh cells and tags', '[2D_OPTIONS]\nDRY_DEPTH 0.001\n', extra_1d=heat,
         files={'tw_heat.heat': heat_cfg})
    deck('tw_heat_nomesh', 'heat mesh overrides without a mesh', '[2D_OPTIONS]\nDRY_DEPTH 0.001\n', extra_1d=heat.replace('tw_heat.heat', 'tw_heat_ok.heat'), mesh=False,
         files={'tw_heat_ok.heat': '[RADIATIVE_FLUXES]\nALBEDO MESH_TAG low 0.3\n'})
    deck('tw_durations', 'duration tokens parse_time_seconds reads as 0 or by its leading hours (U302)',
         '[2D_OPTIONS]\nREPORT_2D_STEP abc\nINFIL_STEP 1h\n[2D_INFILTRATION_OPTIONS]\nINFIL_STEP junk\n')
    deck('tw_output_file', 'OUTPUT_FILE without the HDF5 writer warns (not in metadata mode)', '[2D_OPTIONS]\nOUTPUT_FILE "out dir/r.h5"\n')
    deck('tw_sections_only', 'every 2D section header with no rows', ''.join(f'[{s}]\n' for s in [
        '2D_OPTIONS', '2D_VERTICES', '2D_TRIANGLES', '2D_QUADS', '2D_INITIAL_VELOCITY', '2D_VERTEX_NODE_MAP',
        '2D_TRIANGLE_NODE_MAP', '2D_BOUNDARY_CONDITIONS', '2D_EDGE_CONVEYANCE', '2D_INFILTRATION_OPTIONS',
        '2D_INFILTRATION_DEFAULTS', '2D_INFILTRATION', '2D_MESH_FILE', '2D_INITIAL_QUALITY', '2D_INITIAL_BED_QUALITY',
        '2D_BOUNDARY_QUALITY', 'GW_TRANSPORT_OPTIONS', 'GW_TRANSPORT_PARAMS', 'GW_SORPTION', 'GW_INITIAL_QUALITY',
        'GW_BOUNDARY_QUALITY', 'GW_SOURCES', '2D_COVERAGES', '2D_LOADINGS', '2D_CURB_LENGTH', '2D_SURFACE_OWNERSHIP',
        '2D_AQUIFER_BOUNDARIES', '2D_AQUIFER_OPTIONS', '2D_AQUIFER', '2D_AQUIFER_NODE', '2D_AQUIFER_LINKS']), mesh=False)


if __name__ == '__main__':
    main()
