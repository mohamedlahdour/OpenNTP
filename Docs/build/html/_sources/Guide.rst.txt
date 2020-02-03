User's Guide
============
========================
1. Graphical User Interface
========================

A graphical user interface written in Python programing language has been developed to simplify the use of package OpenRSN.
After starting the software by typing the following command line in a terminal:

    .. code-block:: python

         cd OpenRSN
         $ python3 main.py

A main window (GUI) of the package OpenRSN on an Ubuntu Linux machine will be displayed as in Figure bellow.

.. image:: _images/main.png 

========================
2. Examples & test suite
========================

-----------------------------
2.2. Two-dimensional benchmark
-----------------------------

2.2.1. Creating JSON Input File
-----------------------------

The input data file must be in JSON format. The first method is
to write it directly in **Text editor**  as is shown in the figure bellow, which illustrates an example input file of a simplified nuclear reactor with a quarter  of the x,y square 2D core taken from `[Filho et al.,2002]. <https://doi.org/10.1016/S0168-9274(01)00074-5>`_


    .. code-block:: json
 
        {
          "data": {
            "parameter": {
              "id": 100,
              "Total number of energy groups": 1,
              "Total number of Materials": 4,
              "Total number of X regions": 4,
              "Total number of Y regions": 4,
              "X region thickness per [cm]": [35, 10, 20, 40],
              "Y region thickness per [cm]": [35, 10, 20, 40],
              "Which material goes in each cell":[[4, 4, 4, 4],
                                                  [3, 3, 3, 4],
                                                  [2, 2, 3, 4],
                                                  [1, 2, 3, 4]],
              "XY number of fine meshes in each cell":[[8, 8, 8, 8],
                                                       [8, 8, 8, 8],
                                                       [8, 8, 8, 8],
                                                       [8, 8, 8, 8]],
              "Number of Angular Discretization": 8,
              "The l-order Legendre polonomial": 0,
              "Maximum Iteration": 200,
              "Epsilon Keff": 1.0e-8
            },
            "materials": [
              {
                "id": 1,
                "nom": "a",
                "XSTotal": [0.222589],
                "XSNuFission": [0.00283283],
                "XSScatter Matrix":[[0.220563]],
                "XSChi":  [1.0]
              },
              {
                "id": 2,
                "nom": "b",
                "XSTotal": [0.216566],
                "XSNuFission": [0.0104347],
                "XSScatter Matrix":[[0.210697]],
                "XSChi":  [1.0]
              },
              {
                "id": 3,
                "nom": "c",
                "XSTotal": [0.301439],
                "XSNuFission": [0.000513036],
                "XSScatter Matrix":[[0.296069]],
                "XSChi":  [1.0]
              },
              {
                "id": 4,
                "nom": "d",
                "XSTotal": [0.252250],
                "XSNuFission": [0.0],
                "XSScatter Matrix":[[0.250794]],
                "XSChi":  [0.0]
              }
            ]
          }
        }

2.2.2. Generating JSON Input File
---------------------------------

Another method is to use a set of buttons on the left side of the main
window [ref], these buttons allow users to insert input data automatically without requiring an in-depth knowledge of JSON file syntax. Once the user clicks on **Data Up** button the input JSON file will be automatically generated in the window **Text editor** .

.. image:: _images/insert.png 


2.2.3. Running OpenRSN under a GUI
----------------------------------

The **Run** button is used to running the multi-group scheme, and the figure below shows the values of the multiplication factor as a function of the iteration numbers.

.. image:: _images/runing.png

2.2.4. Geometry Visualization
----------------------------

The **geometry** button allowing to plot in two dimensions the geometry to study. The plotting mode of the geometry is based on the presence of an input file. A depiction of the geometry for the example input file given in sub section `Creating JSON Input File`_ is illustrated in Figure bellow

.. image:: _images/geom.png 

2.2.5. Flux Visualization
-------------------------

The **Plot** button refers to a set of routines programming in fortran and python to plot the scalar flux in space of one or two-dimensional and in each energy group. The figure bellow shows the flux for the example input file given in sub section `Creating JSON Input File`_  with four regions and four materials after clicking on the **Plot** button.

.. image:: _images/flux.png 
2.2.6. Level symmetric gaussian quadrature sets visualization
-------------------------------------------------------------

The level-symmetric quadrature set is used in the Discrete Ordinates (:math:`S_{N}`) method (Lewis and Miller, 1984). The subscript :math:`N`.refers to the number of directions along each axis with half being positive and half negative. The figure below give the weights and angles used for each set in the 1st octant which will be displayed automatically by clicking on the **Ordinate** button

.. image:: _images/ordin.png 

2.2.7. Simple Output
--------------------

The following is the corresponding output to the above case. A brief outline of the output file contents is version and run time information, print of input values of the name list variables, print of relevant parameters after setup, calculation run time parametres method, scalar flux solution and output parameters solution to transport eqaution.
  
    .. code-block:: python

        ********************************************************************************
        ERSN, UNIVERSITY ABDELMALEK ESSAADI FACULTY OF SCIENCES - TETOUAN, MOROCCO
        CODE  DEVELOPED  BY  MOHAMED  LAHDOUR,  PHD  STUDENT
        OpenRSN:         SN  DISCRETE  ORDINATES  METHOD
        DIMENSION:       TWO DIMENSIONS (2D) 
        GEOMETRY:        CARTESIAN
        VERSION NUMBER:  1.2
        VERSION DATE:    20  August  2019
        RAN ON:          2019-08-21 12:45:19.93   (H:M:S)
        ********************************************************************************
                    ----------------------------------------------------------
                              INPUT  PARAMETER - VALUES  FROM  INPUT
                    ----------------------------------------------------------
 
        ENERGY GROUPS NUMBER:                               1
        X REGIONS NUMBER:                                   4
        Y REGIONS NUMBER:                                   4
        MATERIALS NUMBER:                                   4
        SIZE OF EACH X REGION [CM]:          35.00000  10.00000  20.00000  40.00000
        SIZE OF EACH Y REGION [CM]:          35.00000  10.00000  20.00000  40.00000
        NUMBER OF DIRECTION ALONG EACH AXIS:                8
        ORDER LEGENDRE POLYNOMIAL:                          0
        TOTAL NUMBER OF X FINE MESHES:                      8
        TOTAL NUMBER OF Y FINE MESHES:                      8
        CONVERGENCE CRITERION of KEFF AND FLUX:       1.0E-08
 
                    ----------------------------------------------------------
                              CALCULATION  RUN-TIME  PARAMETERS  SN
                    ----------------------------------------------------------
 
        LEVEL  SYMMETRIC  GAUSSIAN  QUADRATURE  SETS: 
 
        N. ORDER          MU             ETA             PSI         WEIGHTS 
 
           1          9.51190E-01     2.18218E-01     2.18218E-01     1.20988E-01
           2          7.86796E-01     2.18218E-01     5.77350E-01     9.07407E-02
           3          7.86796E-01     5.77350E-01     2.18218E-01     9.07407E-02
           4          5.77350E-01     2.18218E-01     7.86796E-01     9.07407E-02
           5          5.77350E-01     5.77350E-01     5.77350E-01     9.25926E-02
           6          5.77350E-01     7.86796E-01     2.18218E-01     9.07407E-02
           7          2.18218E-01     2.18218E-01     9.51190E-01     1.20988E-01
           8          2.18218E-01     5.77350E-01     7.86796E-01     9.07407E-02
           9          2.18218E-01     7.86796E-01     5.77350E-01     9.07407E-02
          10          2.18218E-01     9.51190E-01     2.18218E-01     1.20988E-01
 
        PSEUDO  CROSS  SECTIONS  DATA: 
 
         MATERIAL :   1
 
         GROUP           TOTAL        ABSORPTION      NU*FISSION      SCATTERING      FISSION SPECTRUM
 
           1          2.22589E-01     2.02600E-03     2.83283E-03     2.20563E-01     1.00000E+00
         MATERIAL :   2
 
         GROUP           TOTAL        ABSORPTION      NU*FISSION      SCATTERING      FISSION SPECTRUM
 
           1          2.16566E-01     5.86900E-03     1.04347E-02     2.10697E-01     1.00000E+00
         MATERIAL :   3
 
         GROUP           TOTAL        ABSORPTION      NU*FISSION      SCATTERING      FISSION SPECTRUM
 
           1          3.01439E-01     5.37000E-03     5.13036E-04     2.96069E-01     1.00000E+00
         MATERIAL :   4
 
         GROUP           TOTAL        ABSORPTION      NU*FISSION      SCATTERING      FISSION SPECTRUM
 
           1          2.52250E-01     1.45600E-03     0.00000E+00     2.50794E-01     0.00000E+00
 
                    ----------------------------------------------------------
                                NORMALIZED SCALAR  FLUX  SOLUTION
                    ----------------------------------------------------------
 
        FLUXES  PER  MESH  PER  ENERGY  GROUP:
 
        M E S H       G R O U P 1
 
                     1            2            3            4            5            6            7            8
            1   1.00000E+00  8.99018E-01  7.96162E-01  6.93054E-01  4.64773E-01  2.39120E-01  1.17909E-01  3.51124E-02
            2   9.01165E-01  7.97144E-01  6.94076E-01  5.98010E-01  3.94987E-01  2.00672E-01  1.00336E-01  3.05219E-02
            3   7.99305E-01  6.95147E-01  5.78156E-01  4.86665E-01  3.19286E-01  1.62066E-01  8.34442E-02  2.61713E-02
            4   6.95459E-01  5.98931E-01  4.86670E-01  4.05703E-01  2.69057E-01  1.41666E-01  7.46702E-02  2.40013E-02
            5   4.66602E-01  3.95815E-01  3.19399E-01  2.69137E-01  1.90544E-01  1.12075E-01  6.16261E-02  2.09120E-02
            6   2.40387E-01  2.01169E-01  1.62285E-01  1.41879E-01  1.12114E-01  7.75854E-02  4.62275E-02  1.64931E-02
            7   1.18727E-01  1.00733E-01  8.35693E-02  7.48131E-02  6.17163E-02  4.62580E-02  2.92000E-02  1.04020E-02
            8   3.53836E-02  3.06810E-02  2.62660E-02  2.40611E-02  2.09463E-02  1.65158E-02  1.04078E-02  3.69251E-03
 
                    ----------------------------------------------------------
                      OUTPUT  PARAMETER - SOLUTION  TO  TRANSPORT  EQUATION
                    ----------------------------------------------------------
 
        K-EFF                    =      0.968367
        N. OUTER ITERATIONS      =          7331
        TOTAL INNER ITERATIONS   =            15
        TOTAL EXECUTION TIME     =    0:00:01.65   (H:M:S)
 
        ********************************************************************************

-----------------------------
2.3. Slab
-----------------------------
Setting up input file for slab geometry in two energy groups with isotropic scattering source.

    .. code-block:: json

        { 
          "data": { 
            "parameter": { 
              "id": 100,
              "Total number of energy groups": 2,
              "Total number of Materials": 2,
              "Total number of regions": 3,
              "Which material goes in each region": [2, 1, 2],
              "Size for each material per [cm]": [5.630757, 9.726784, 5.630757],
              "Number of fine meshes": [50, 50, 50],
              "Number of Angular Discretization": 8,
              "The l-order Legendre polonomial": 0,
              "Maximum Iteration": 200,
              "Epsilon Keff": 1e-8
            }, 
            "materials": [
              { 
                "id": 1, 
                "nom": "material 1",
                "XSTotal": [0.88721, 2.9727],
                "XSNuFission": [0.00209, 0.07391],
                "XSScatter Matrix":[[[0.838920, 0.04635],
                                     [0.000767, 2.91830]]],
                "XSChi":  [1.0, 0.0]
              },
              { 
                "id": 2, 
                "nom": "material 2",
                "XSTotal": [0.88798, 2.9865],
                "XSNuFission": [0, 0],
                "XSScatter Matrix":[[[0.83975, 0.04749],
                                     [0.000336, 2.96760]]],
                "XSChi":  [0, 0]
              }
            ]  
          }  
        }

Geometry in a one-dimensional slab

.. image:: _images/SlabG.png 

Flux in a one-dimensional slab

.. image:: _images/SlabF.png 

-----------------------------
2.4. Cylinder
-----------------------------

An example for cylindrical  infinite  cell equivalent to the **TRIGA MARK-II** research reactor pin cell is presented here by using 7 energy groups

    .. code-block:: json

        { 
          "data": { 
            "parameter": { 
              "id": 100,
              "Total number of energy groups": 7,
              "Total number of Materials": 5,
              "Total number of regions": 5,
              "Which material goes in each region": [1, 2, 3, 4, 5],
              "Ray for each region per [cm]": [0.3175, 1.82769, 1.83150, 1.88230, 2.285814144],
              "Number of fine meshes": [5, 5, 5, 5, 5],
              "Number of Angular Discretization": 0,
              "The l-order Legendre polonomial": 0,
              "Maximum Iteration": 200,
              "Epsilon Keff": 1.0e-8
            }, 
            "materials": [
              { 
                "id": 1, 
                "nom": "material 1",
                "XSTotal": [0.297551431, 0.288909664, 0.290306468, 0.286637159, 0.295583239, 0.326837471, 0.155639234],
                "XSNuFission": [0, 0, 0, 0, 0, 0, 0],
                "XSScatter Matrix":[[[0.269680893, 0.021221626, 0, 0, 0, 0, 0], 
                                     [0.015676686, 0.263059691, 0.005057992, 0, 0, 0, 0], 
                                     [0, 0.013507066, 0.273815991, 0.001186857, 0, 0, 0],
                                     [0, 0, 0.014034288, 0.270365363, 0.002092507, 0, 0], 
                                     [0, 0, 0, 0.00079017, 0.292969338, 0,  0], 
                                     [0, 0, 0, 0, 0.001078437, 0.32539832 , 0], 
                                     [0, 0, 0, 0, 0.000007539, 0.011162732, 0.14376993]]],
                "XSChi":  [0, 0, 0, 0, 0, 0, 0]
              },
              { 
                "id": 2, 
                "nom": "material 2",
                "XSTotal": [2.994544,  1.408450959, 0.817921932, 0.633026483, 0.650866502, 0.515434679, 0.243844582],
                "XSNuFission": [0.362087986, 0.179080763, 0.09568437, 0.03889374, 0.013653221, 0.001011529, 0.00163296],
                "XSScatter Matrix":[[[2.739733851, 0.036287563, 0.014790886, 0, 0, 0, 0],
                                           [0.038125613, 1.255653823, 0.014136713, 0, 0, 0, 0], 
                                           [0.166974328, 0.274556127, 0.320581789, 0.00109001, 0, 0, 0], 
                                           [0.036264407, 0.11429219, 0.651920757, -0.19634513, 0.001532029, 0, 0], 
                                           [0.003270334, 0.00846593, 0.051970962, 0.052892578, 0.512276254, 0, 0],
                                           [0.000000452, 0.000001176, 0.000010132, 0.000008685, 0.162129535, 0.352055605, 0], 
                                           [0, 0, 0.000000177, 0.000000059, 0.001637495, 0.098331459, 0.142495556]]],
                "XSChi":  [0, 0, 0, 0, 0.000390, 0.126214, 0.872105]
              },
              { 
                "id": 3, 
                "nom": "material 3",
                "XSTotal": [0, 0, 0, 0, 0, 0, 0],
                "XSNuFission": [0, 0, 0, 0, 0, 0, 0],
                "XSScatter Matrix":[[[0, 0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 0, 0, 0, 0], 
                                     [0, 0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 0, 0, 0, 0]]],
                "XSChi":  [0, 0, 0, 0, 0, 0, 0]
              },
              { 
                "id": 4, 
                "nom": "material 4",
                "XSTotal": [1.153239138, 1.012997908, 0.947913799, 0.908546838, 0.966825787, 0.46508365, 0.201018163],
                "XSNuFission": [0, 0, 0, 0, 0, 0, 0],
                "XSScatter Matrix":[[[0.810739629, 0.066072482, 0, 0, 0, 0, 0], 
                                     [0.070201393, 0.778857964, 0.018049333, 0, 0, 0, 0], 
                                     [0.000005432, 0.055959691, 0.806723946, 0.003927282, 0, 0, 0],
                                     [0, 0, 0.061596098, 0.790102323, 0.005981459, 0, 0],
                                     [0, 0, 0, 0.003658891, 0.950937049, 0 , 0], 
                                     [0, 0, 0, 0, 0.005035727, 0.459136241, 0], 
                                     [0, 0, 0, 0, 0.000075732, 0.014066679, 0.185277698]]],
                "XSChi":  [0, 0, 0, 0, 0, 0, 0]
              },
              { 
                "id": 5, 
                "nom": "material 5",
                "XSTotal": [3.102999797, 1.73516897, 0.980823407, 0.642418097, 0.567747919, 0.419879925, 0.201975672],
                "XSNuFission": [0, 0, 0, 0, 0, 0, 0],
                "XSScatter Matrix":[[[2.631785406, 0.437581882, 0.00907283, 0, 0, 0, 0],
                                     [0.609582022, 1.041127377, 0.071965385, 0, 0, 0, 0],
                                     [0.141504305, 0.538385484, 0.292249692, 0.002089256, 0, 0, 0], 
                                     [0.036338638, 0.123430676, 0.814518036, -0.336757332, 0.001951418, 0, 0], 
                                     [0.003425563, 0.009933065, 0.061937554, 0.066020159, 0.425993565, 0, 0], 
                                     [0.000000363, 0.000001996, 0.000009253, 0.00000889,  0.203705679, 0.216356046, 0], 
                                     [0, 0, 0.000000253, 0, 0.001995558, 0.110723501, 0.089154896]]],
                "XSChi":  [0, 0, 0, 0, 0, 0, 0]
              }
            ]  
          }  
        }

Geometry in a two-dimensional TRIGA Reactor 

.. image:: _images/TrigaG.png 

The infinite cell in OpenMC `OpenMC <https://openmc.readthedocs.io/en/stable/>`_ is represented by hexagonal cell with reflective boundaries. The infinite multiplication factor values ​​obtained in `OpenRSN <https://openrsn.readthedocs.io/en/latest/index.html>`_ and `OpenMC <https://openmc.readthedocs.io/en/stable/>`_ are shown in Table below.

.. table:: Calculate infinite multiplication factor :math:`k_{inf}`.

    +----------------------+---------------+
    | Surface              |:math:`k_{inf}`|
    +======================+===============+
    | OpenRSN              | 1.403180      |
    |                      |               |
    +----------------------+---------------+
    |                      |               |
    | OoenMC               | 1.417226      |
    +----------------------+---------------+

-----------------------------
2.5. Sphere
-----------------------------

-----------------------------
2.6. Pin Cell
-----------------------------