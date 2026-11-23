This scatter1D_readme.txt file was generated on 20200814 by Bryan Barnes

-----------------------------------------------------------
GENERAL INFORMATION
-----------------------------------------------------------

1. Title of Dataset 
    Optical scattering measurements and simulation data for one dimensional 
	(1-D) patterned periodic sub-wavelength features 

2. Format of the files
    Comma-separated values

3. Author Information

  Principal Investigator Contact Information
   Name:        Bryan M. Barnes
   Institution: Nanoscale Device Characterization Division,
                Physical Measurement Laboratory,
                National Institute of Standards and Technology (NIST)
   Address:     100 Bureau Dr MS 8423, Gaithersburg, MD 20878
   Email:       bmbarnes@nist.gov
   
  Associate Contact Information 
   Name:        Mark-Alexander Henn
   Institution: Applied and Computational Mathematics Division,
                Information Technology Laboratory,
                National Institute of Standards and Technology (NIST)
   Address:     100 Bureau Dr, Gaithersburg, MD 20878
   Email:       mark.henn@nist.gov 

  Alternate Contact Information
   Name:        Richard M. Silver
   Institution: Nanoscale Device Characterization Division,
                Physical Measurement Laboratory,
                National Institute of Standards and Technology (NIST)
   Address:     100 Bureau Dr MS 8423, Gaithersburg, MD 20878
   Email:       silver@nist.gov

4. Date of data collection 
    20080714-20080805 (experimental)     
    20080423-20081031 (simulation)

5. Geographic location of data collection (if relevant. include LAT/LONG 
   coordinates, if relevant):
    not relevant

6. Information about any NON-NIST funding sources that supported the collection of 
   the data:
   
-----------------------------------------------------------
SHARING/ACCESS INFORMATION
-----------------------------------------------------------
7.  Other than the NIST statements for Copyright, Fair Use, and Licensing found 
    at
https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software, 
    list the other license or restrictions which are placed on this data:
     none

8. Are there any restrictions or guidelines on how to use the data, 
   (ex. information about access restrictions based on privacy, security,
   or other policies)?
     none

9. Is there a documentary standard that applies to this dataset? 
     none

10. Links to publications that cite or use the data:
     https://doi.org/10.1117/12.816569  
     https://doi.org/10.1117/12.827676
     https://doi.org/10.1364/AO.51.006196
     https://doi.org/10.1088/1361-6501/aa5586
     https://doi.org/10.1117/12.2551504

11. Links to other publicly accessible locations of the data:
     none

12. Links/relationships to ancillary data sets:
     none

13. Was data derived from another source? (example, Open FEMA)
     no

14. Recommended citation for the data 
     B.M. Barnes, M.-A. Henn, R.M. Silver, (2020), Optical Scattering 
     Measurements and Simulation Data for 1-D Patterned Periodic Sub-
     Wavelength Features, National Institute of Standards and Technology, 
     https://doi.org/10.18434/mds2-2290 (Accessed MONTH DAY, YEAR).


-----------------------------------------------------------
DATA & FILE OVERVIEW
-----------------------------------------------------------

15. File List
   A. Filename:             Pi_L100P300_sim_All.csv        
      Short description:    2566 x 3 matrix of parametric values used as inputs 
                            the scattering code for simulation.  Values are in 
                            nanometers.
      
   B. Filename:             Phi_L100P300_sim_All.csv
      Short description:    2566 x 84 matrix of simulated intensities from the
                            scattering off a periodic structure parametrized using 
                            a double-trapezoid cross-section.  These intensities 
                            are normalized and unitless as described below.
        
   C. Filename:             Imean_L100P300_exp.csv
      Short description:    9 x 84 matrix of measurements determined from repeated 
                            measurements of the intensities scattered off the nine 
                            measured "L100P300" targets, specifically within nine 
                            dies on a single wafer.  Each row of intensities 
                            corresponds to a single die.  These intensities are 
                            normalized and unit-less.
                            
   D. Filename:             Istd_L100P300_exp.csv
      Short description:    9 x 84 matrix of Type A measurement errors propagated 
                            in the measurement and normalization of 
                            "Imean_L100P300_exp.csv" Each row of errors 
                            corresponds to an individual die.

   E. Filename:             labels_L100P300.csv
      Short description:    3 x 84 matrix of labels assignable to the 84 columns.
                            Rows correspond to linear polarization state, plane of 
                            incidence, and angle of incidence (in degrees).
                                                        
   F. Filename:             y_L100P300_Die1_exp.csv        
      Short description:    1 x 84 vector formatted for use with our group's 
                            published regression software, M.o.R. 
                            http://doi.org/10.18434/T4/1502429 
                            It is a subset of "Imean_L100P300_exp.csv"

   G. Filename:             invV_L100P300_Die1_exp.csv
      Short description:    84 x 84 matrix that is formatted for use with M.o.R.  
                            It is the inverse of the V matrix in the regression.  
                            It corresponds to the same die as that for 
                            "y_L100P300_Die1.csv".  Unitless.  
            
   H. Filename:             Pi_L100P300_sim_467.csv        
      Short description:    467 x 3 matrix of parametric values used as inputs 
                            the scattering code for simulation.  Values in 
                            nanometers. Formatted for use in M.o.R., this is a 
                            subset of the larger dataset.
      
   I. Filename:             Phi_L100P300_sim_467.csv
      Short description:    467 x 84 matrix of simulated intensities from the
                            scattering off a periodic structure parameterized using 
                            a double-trapezoid cross-section.  Normalized and 
                            unitless.  Formatted for use in M.o.R., this is a 
                            subset of the larger dataset.
							
   J. Filename:             Pi_L100P300_sim_140.csv        
      Short description:    140 x 3 matrix of parametric values used as inputs the 
                            scattering code for simulation.  Values in nanometers.
                            This is a subset of the larger dataset formatted for use
                            in M.o.R., sized for use in machine learning.
      
   K. Filename:             Phi_L100P300_sim_140.csv
      Short description:    140 x 84 matrix of simulated intensities from the
                            scattering off a periodic structure parametrized using 
                            a double-trapezoid cross-section.  Normalized and 
                            unitless. This is a subset of the larger dataset
                            formatted for use in M.o.R., sized for use in machine 
                            learning.
	
   L. Filename:             Pi_L100P300_sim_Comparison.png
      Short description:    Image illustrating the sampling of the full three-
                            parameter space using the full set of 2566 simulations                     
                            and the smaller set of 467 simulations.
							
   K. Filename:             Pi_L100P300_sim_20_through_140.png
      Short description:    Image illustrating the 140 x 3 sampling of the full three-
                            parameter space as it developed by multiples of 20 using 
                            a Halton sequence.
							
16. Relationship between files:        
     --  Filenames with the suffix "exp" are experimental data, 
         suffix "sim" are simulation data.
     --  Pi_L100P300_sim_All.csv are input values used to determine 
         Phi_L100P300_sim_All.csv from electromagnetic simulations.
     --  labels_L100P300.csv describe the polarization axis, incidence plane, and 
         angle of incidence for the following files: 
          * Phi_L100P300_sim_All.csv
          * Phi_L100P300_sim_467.csv
          * Phi_L100P300_sim_140.csv
          * Imean_L100P300_exp.csv
          * Istd_L100P300_exp.csv
          * y_L100P300_Die1_exp.csv
     --  A subset of data is provided as separate files for direct use in our  
          group's "M.o.R." radial basis function fitting program 
          http://doi.org/10.18434/T4/1502429 :
           *  Both (Pi_L100P300_sim_100.csv and Phi_L100P300_sim_100.csv) and
	      (Pi_L100P300_sim_467.csv and Phi_L100P300_sim_467.csv) are subsets of
              Pi_L100P300_sim_All.csv and Phi_L100P300_sim_All.csv with less 
              sampling density with matrix formats ready for use in M.o.R.
	      A visual comparison between the 467 vs full set is provided as 
	      "Pi_L100P300_sim_Comparison.png".
           *  y_L100P300_Die1.csv is the first row of Imean_L100P300_exp.csv .  
           *  invV_L100P300_Die1.csv uses the first row of uncertainties in 
              y_L100P300_Die1.csv to estimate variances which are then the 
	      diagonal of the 84 x 84 matrix.
    
17. Additional related data collected that was not included in the current data 
    package:
     Raw data from each L100P300 target and from silicon wafers were collected
     under illumination.  Raw dark count data were also collected in the 
     absence of illumination.  Glare data were gathered with the illumination
     on in the absence of a sample to capture excess light reaching the 
     camera reflected off the objective lens.  
     
18. Are there multiple versions of the dataset? yes/no
     no

-----------------------------------------------------------
METHODOLOGICAL INFORMATION
-----------------------------------------------------------

19. Description of methods used for collection/generation of data: 
     These experimental data were collected using Angle-resolved 
     scatterfield microscopy. In this configuration, illumination at the sample is 
     governed by manipulation of Kohler illumination at the conjugate to the back 
     focal plane (CBFP) of the objective lens.  A circular aperture in the CBFP 
     determines both the angular resolution of light at the sample (due to 
     diameter) and also both the incident plane and angle of illumination at the 
     sample plane (due to aperture position).  Expressed as numerical aperture 
     NA, the angular resolution is 0.13 NA.  A linear polarizer very close to 
     the CBFP governs the polarization state of this narrowed illumination beam.
         The target itself is utilized to define the coordinates used in the
     experiment:
     
      Y                             | | | | | | |    
      ^                             | | | | | | |     (schematic of 
      |       (coordinates,         | | | | | | |      periodic structure,
      +-->X    Z,out of page)       | | | | | | |      top view)
     
     Four types of angle scans have been performed using 21 angles each.  Two 
     of the scans have the polarization pass axis aligned with the Y direction
     above, two aligned with the X direction.  Likewise, two of the scans have the
     plane of incidence across the lines (i.e. in the X-Z plane) and two have the
     plane of incidence along the lines (i.e. in the Y-Z plane).  These are mixed 
     such that each of the four scans has a different combination of polarization
     pass axis and plane of incidence.  Data from these scans can be concatenated
     to form the vectors of length 84 that are presented here.
       Additional data are collected contemporaneously on bare silicon wafers.  
     In addition, with the illumination off, the dark current of the charge-
     coupled device (CCD) camera is measured.  Also, with the sample removed and 
     the illumination on, the intensity of angle-resolved glare is recorded at the
     camera. This glare is light that does not reach the sample but rather 
     reflects off the physical optics and into the collection path.  This glare 
     intensity by definition also includes the effects of the dark current.
        For the simulations, three parametric values were used as floating parameters
     that defined the widths of stacked double trapezoids.  These trapezoids of equal 
     height are the geometric inputs for a rigorous coupled-wave analysis (RCWA) 
     of the optical scattering off the nominal 100 nm linewidth structures.  
     Specifically, input files fixed the total height of the pair of trapezoids to 
     229 nm with a fixed pitch of 301.5 nm.  The optical constants for simulation 
     were fixed to (4.6959 + 0.0802i) for both the line and substrate.  
   
20. Methods for processing the data: 
     The raw intensity data collected from the scattering is biased by the 
     ever-present dark current and glare. By subtracting the measured glare from 
     both the bare silicon wafer data and the measured scattering from the
     periodic structure and taking their ratio, one learns the relative strength
     of the scattering.  Furthermore, as the reflectivity of silicon 
     is easily calculated using Fresnel equations, the measured intensity 
     can be normalized relative to the incident intensity I_0 = 1, an assumption 
     underlying reflectivity and also our electromagnetic modeling.  
        A full description of these experimental methods is provided in 
     https://doi.org/10.1088/0957-0233/22/2/024003      

21. Instrument- or software-specific information needed to interpret the data:
     Illumination NA = 0.13, collection NA = 0.95
	 These are the settings both in experiment and simulation.
	 The illumination wavelength was 450 nm.

22. Standards and calibration information, if appropriate:
     none

23. Environmental/experimental conditions:
     Samples were measured within a temperature controlled, clean room 
	 environment.

24. Describe any quality-assurance procedures performed on the data:
     none

25. People involved with sample collection, processing, analysis and/or 
    submission:
     * Rick Silver supervised both data analysis and design of experiment.
     * Bryan Barnes performed data collection, processing, additional analysis, 
       and submission.
     * Mark-Alexander Henn performed additional analysis.
     * The rigorous couple-wave analysis (RCWA) code and simulated imaging code
       were developed in-house by Thomas A. Germer of NIST.


-----------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: Pi_L100P300_sim_All.csv
-----------------------------------------------------------
   
26. Number of variables:
     3 

27. Number of cases/rows: 
     2566

28. Variable List
    A. Name: widthBottom
       Description: parametric value for the bottom width of a double trapezoid 
	   used as inputs a rigorous coupled-wave analysis (RCWA) scattering code 
	   for simulating the scattering off the nominal 100 nm linewidth 
	   structures with nominal 300 nm pitch.  Units are nanometers.  

    B. Name: widthMiddle
       Description: parametric value for the middle width of a double trapezoid
	   used as inputs to the RCWA scattering code.  Units are nanometers.
       
    C. Name: widthTop
       Description: parametric value for the top width of a double trapezoid 
	   used as inputs to the RCWA scattering code.  Units are nanometers.

29. Missing data codes:
     none
     
30. Specialized formats of other abbreviations used
     none

-----------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: Phi_L100P300_sim_All.csv
-----------------------------------------------------------

26. Number of variables:
     84 

27. Number of cases/rows: 
     2566

28. Variable List
    A. Names:       intensitySimPolXAngXnn for nn = 01,...,21
       Description: Normalized unitless intensities for simulated, linearly 
        polarized illumination across the features (X pol.).  The incident 
        plane is aligned across the features (X axis).

    B. Name:        intensitySimPolYAngXnn for nn = 22,...,42
       Description: Normalized unitless intensities for simulated, linearly 
        polarized illumination along the features (Y pol.).  The incident 
        plane is aligned across the features (X axis).

    C. Name:        intensitySimPolXAngYnn for nn = 43,...,63
       Description: Normalized unitless intensities for simulated, linearly 
        polarized illumination across the features (X pol.).  The incident 
        plane is aligned along the features (Y axis).

    D. Name:        intensitySimPolYAngYnn for nn = 64,...,84
       Description: Normalized unitless intensities for simulated, linearly 
        polarized illumination along the features (Y pol.).  The incident 
        plane is aligned along the features (Y axis).

29. Missing data codes:
     none

30. Specialized formats of other abbreviations used
     none
     
-----------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: Imean_L100P300_exp.csv
-----------------------------------------------------------

26. Number of variables:
     84

27. Number of cases/rows: 
     9

28. Variable List
    A. Names:       intensityExpPolXAngXnn for nn = 01,...,21
       Description: Normalized unitless intensities for measured linearly 
        polarized illumination across the features (X pol.).  The incident 
        plane is aligned across the features (X axis).

    B. Name:        intensityExpPolYAngXnn for nn = 22,...,42
       Description: Normalized unitless intensities for measured, linearly 
        polarized illumination along the features (Y pol.).  The incident 
        plane is aligned across the features (X axis).

    C. Name:        intensityExpPolXAngYnn for nn = 43,...,63
       Description: Normalized unitless intensities for measured, linearly 
        polarized illumination across the features (X pol.).  The incident 
        plane is aligned along the features (Y axis).

    D. Name:        intensityExpPolYAngYnn for nn = 64,...,84
       Description: Normalized unitless intensities for measured, linearly 
        polarized illumination along the features (Y pol.).  The incident 
        plane is aligned along the features (Y axis).


      Rows 1-9 correspond to the following targets as identified in a 
       previous publication (https://doi.org/10.1117/12.777131) and summarized 
       here:
       +========+========+========+
       |  Die 1 |  Die 2 |  Die 3 |
       | (1,-1) | (1,0)  |  (1,1) |
       +========+========+========+ 
       |  Die 4 |  Die 5 |  Die 6 |
       | (0,-1) | (0,0)  |  (0,1) |
       +========+========+========+
       |  Die 7 |  Die 8 |  Die 9 |
       |(-1,-1) | (-1,0) | (-1,1) |
       +========+========+========+
     There are physical differences among these die as they were pattered 
     using a focus/exposure matrix.

29. Missing data codes:
     none
	 
30. Specialized formats of other abbreviations used
     none
     
-----------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: Istd_L100P300_exp.csv
-----------------------------------------------------------

26. Number of variables:
     84

27. Number of cases/rows: 
     9

28. Variable List
    A. Names:       uncertExpPolXAngXnn for nn = 01,...,21
       Description: Type A propagated uncertainties (coverage factor = 1) 
        based upon measurements of bare silicon, the L100P300 target, and  
        the dark current with multiple repetitions. Illumination is 
        linearly polarized across the features (X pol.).  The incident 
        plane is aligned across the features (X axis).
		
    B. Names:       uncertExpPolYAngXnn for nn = 22,...,42             
       Description: Type A propagated uncertainties (coverage factor = 1) 
        based upon measurements of bare silicon, the L100P300 target, and  
        the dark current with multiple repetitions. Illumination is 
        linearly polarized along the features (Y pol.).  The incident 
        plane is aligned across the features (X axis).

    C. Names:       uncertExpPolXAngYnn for nn = 43,...,63             
       Description: Type A propagated uncertainties (coverage factor = 1) 
        based upon measurements of bare silicon, the L100P300 target, and  
        the dark current with multiple repetitions. Illumination is 
        linearly polarized across the features (X pol.).  The incident 
        plane is aligned along the features (Y axis).

    C. Names:       uncertExpPolYAngnn for nn = 64,...,84             
       Description: Type A propagated uncertainties (coverage factor = 1) 
        based upon measurements of bare silicon, the L100P300 target, and  
        the dark current with multiple repetitions. Illumination is 
        linearly polarized along the features (Y pol.).  The incident 
        plane is aligned along the features (Y axis).

29. Missing data codes:
     none
	 
30. Specialized formats of other abbreviations used
     none
     
-----------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: labels_L100P300.csv
-----------------------------------------------------------

26. Number of variables:
     3

27. Number of cases/rows: 
     84

28. Variable List
    A. Name:        planePolarization
	Description: either "PolAxX" or "PolAxY" corresponding to the two orthogonal 
        linear polarization states X polarization and Y polarization,
        respectively.

    B. Name:        planeIncidence
       Description: either "AngleX" or "AngleY" corresponding to the two 
	   orthogonal planes defined above, X-Z or Y-Z, respectively.
       
    C. Name:        angleIncidence
       Description: Numerical value in degrees of the angle between the 
                    illumination and the surface normal.  The plane of 
                    incidence is defined in the second column.  Note, values 
                    repeat but the experimental conditions change as shown in
                    columns 1 and 2 to yield 84 measurement conditions.

     The 84 measurement conditions (i.e. rows) are described as variables 
     (columns) in Phi_L100P300_sim_All.csv above.


29. Missing data codes:
     none 
	 
30. Specialized formats of other abbreviations used
     none
     
-----------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: y_L100P300_Die1_exp.csv
-----------------------------------------------------------

26. Number of variables:
     1 

27. Number of cases/rows: 
     84

28. Variable List
    A. Name:        y
       Description: This file is the first row of Imean_L100P300_exp.csv
       separated out for use with our group's "M.o.R." radial basis 
       function fitting program.  
       The 84 measurement conditions (i.e. rows) are described as 
       variables (columns) in Imean_L100P300_exp.csv above.
            
29. Missing data codes:
     none
	 
30. Specialized formats of other abbreviations used
     none

-----------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: invV_L100P300_Die1_exp.csv
-----------------------------------------------------------

26. Number of variables: 
     1 (see Variable List)

27. Number of cases/rows: 
     84

28. Variable List
    A. Name:        invV
       Description: This variable populates the diagonal term of this 
        84 x 84 matrix, with all other entries in the matrix equal to 0.
        Variances are based on the first row of Imean_L100P300_exp.csv
        and formatted for use with our group's "M.o.R." radial basis 
        function fitting program. 

29. Missing data codes:
     none

30. Specialized formats of other abbreviations used
     none
     
-----------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: Pi_L100P300_sim_467.csv
-----------------------------------------------------------

26. Number of variables:
     3

27. Number of cases/rows: 
     467

28. Variable List
     A. Name: widthBottom
       Description: parametric value for the bottom width of a double trapezoid 
        used as inputs a rigorous coupled-wave analysis (RCWA) scattering code 
	for simulating the scattering off the nominal 100 nm linewidth 
	structures with nominal 300 nm pitch.  Units are nanometers.  

    B. Name: widthMiddle
       Description: parametric value for the middle width of a double trapezoid
        used as inputs to the RCWA scattering code.  Units are nanometers.
       
    C. Name: widthTop
       Description: parametric value for the top width of a double trapezoid 
        used as inputs to the RCWA scattering code.  Units are nanometers.
	  
	  Rows have been removed from Pi_L100P300_sim_All.csv to 
      facilitate use with our group's "M.o.R." radial basis 
      function fitting program, which cannot easily process 
      the full number (2566) of matrix entries. 
	 	 
29. Missing data codes:
     none
	 
30. Specialized formats of other abbreviations used
     none
     
-----------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: Phi_L100P300_sim_467.csv
-----------------------------------------------------------

26. Number of variables:
     84

27. Number of cases/rows: 
     467

28. Variable List
     see "Variable List for Phi_L100P300_sim_All.csv" above
      Rows have been removed from Phi_L100P300_sim_All.csv to 
      facilitate use with our group's "M.o.R." radial basis 
      function fitting program, which cannot easily process 
      the full number (2566) of matrix entries. 

29. Missing data codes:
     none
	 
30. Specialized formats of other abbreviations used
     none
	 
-----------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: Pi_L100P300_sim_467.csv
-----------------------------------------------------------

26. Number of variables:
     3

27. Number of cases/rows: 
     140

28. Variable List
     A. Name: widthBottom
       Description: parametric value for the bottom width of a double trapezoid 
        used as inputs a rigorous coupled-wave analysis (RCWA) scattering code 
	for simulating the scattering off the nominal 100 nm linewidth 
	structures with nominal 300 nm pitch.  Units are nanometers.  

    B. Name: widthMiddle
       Description: parametric value for the middle width of a double trapezoid
        used as inputs to the RCWA scattering code.  Units are nanometers.
       
    C. Name: widthTop
       Description: parametric value for the top width of a double trapezoid 
        used as inputs to the RCWA scattering code.  Units are nanometers.
	  
       Rows have been removed from Pi_L100P300_sim_All.csv to facilitate use 
       not only with our group's "M.o.R." radial basis function fitting program,
       but also to enable the use of machine learning techniques such as 
       Gaussian process regression.
	 	 
29. Missing data codes:
     none
	 
30. Specialized formats of other abbreviations used
     none
     
-----------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: Phi_L100P300_sim_100.csv
-----------------------------------------------------------

26. Number of variables:
     84

27. Number of cases/rows: 
     467

28. Variable List
     see "Variable List for Phi_L100P300_sim_All.csv" above
       Rows have been removed from Pi_L100P300_sim_All.csv to facilitate use 
       not only with our group's "M.o.R." radial basis function fitting program,
       but also to enable the use of machine learning techniques such as 
       Gaussian process regression.

29. Missing data codes:
     none
	 
30. Specialized formats of other abbreviations used
     none
