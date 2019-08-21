.. OpenRSN documentation master file, created by
   sphinx-quickstart on Sun Aug 18 20:09:06 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation: OpenRSN
===================================

Neutron Transport Package with Graphical User Interface

Neutron Transport Package **OpenRSN** (Open Neutron Transport Package from the Radiations and Nuclear Systems Group), is an open-source code written in FORTRAN90 for a pedagogical purpose to solve the steady-state multigroup neutron transport equation using either:

- Collision Probablity Method (CP) in One-Dimensional for Slab, Cylindrical or Spherical geoemtry. 
- Discrete Ordinate Method (:math:`S_{N}`) in One or Two-Dimensional for Cartesian Geometry.
- Method of Characteristics (MOC) in One-Dimensional for Slab Geometry.

The code, including the graphical user interface is developed and maintained by `Mohamed LAHDOUR <https://github.com/mohamedlahdour>`_ (PhD student) and Prof. `Tarek EL BARDOUNI  <https://github.com/tarekbardouni>`_ from `University Abdelmalek Essaadi Tetouan Morocco <http://www.fst.ac.ma/site/>`_ .

The code with its GUI is developed and maintained by `Mohamed LAHDOUR <https://github.com/mohamedlahdour>`_ (PhD student) and Prof. `Tarek EL BARDOUNI  <https://github.com/tarekbardouni>`_ from `University Abdelmalek Essaadi Tetouan Morocco <http://www.fst.ac.ma/site/>`_ .

**OpenRSN’s main features are:**

* free & open source software with a pedagogical purpose.
* solve the steady-state multigroup neutron transport equation in one, or two spatial dimensions.
* solve the steady-state multigroup neutron transport equation  in a multiplicative medium with isotropic and anisotropic scatternig source.
* simple framework to add and test new algorithms.
* provided with a graphical user interface written in Python programing language which has been developed to simplify the use of **OpenRSN**.
* computed results: :math:`k_{eff}` and fluxes.



.. admonition:: Recommended publication for citing
   :class: tip

   Lahdour, M., Bardouni, T. E., Chakir, E., Benaalilou, K., Mohammed, M., Bougueniz, H.,
   and Yaakoubi, H. E., "`Ntp-ersn: A new package for solving the multigroup neutron
   transport equation in a slab geometry. <https://doi.org/10.1016/j.apradiso.2018.12.004>`_,"
   *Applied Radiation and Isotopes*, 145:73 – **84** (2019a).


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Installation
   Theory
   Guide
   License
   Help
   Publications
   developers
   



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
