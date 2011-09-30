# Simple counting experiment, with one signal and one background process
#imax 1  number of channels
#jmax 1  number of backgrounds 
#kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------
shapes ggH CMS_hwwlvqq_ele hwwlvqq_ele.input.root  w:signal
shapes VBF CMS_hwwlvqq_ele hwwlvqq_ele.input.root  w:signal 
shapes background CMS_hwwlvqq_ele hwwlvqq_ele.input.root w:background
###shapes data_obs CMS_hwwlvqq_ele hwwlvqq_ele.input.root w:data_obs
shapes data_obs CMS_hwwlvqq_ele hwwlvqq_ele.input.root w:dataset_obs
------------
bin         CMS_hwwlvqq_ele
observation <dummyobs>
------------
bin                CMS_hwwlvqq_ele	  CMS_hwwlvqq_ele	          CMS_hwwlvqq_ele
process       ggH       		  				  VBF                 background
process         -1                                                          0                        1        
rate	    <dummy1>
------------
lumi		lnN	1.045			1.045			1.0
pdf_ggH	         <dummypdfggH>
pdf_qqH	       <dummypdfqqH>
QCDscale_ggH				 <dummyggH>
QCDscale_qqH			 <dummyVBF>
#theory_gamma                      <dummygammaBW>
CMS_trigger_ele	lnN	1.02	1.02	1.0	
CMS_eff_ele	lnN	1.008	1.008	1.0
CMS_scale_ele	lnN	1.01	1.01	1.0
#CMS_reco_ele	lnN	1.015	1.015	1.0
CMS_scale_j	<dummyJES>
CMS_eff_b 	<dummybeff>
CMS_hwwlvqq_pu	lnN		1.02	        1.02		1.0			      	
CMS_hwwlvqq_qgsep   lnN	1.046            1.046          1.0
CMS_hwwlvqq_bkgelep0   <dummybnorm>
CMS_hwwlvqq_bkgelep1   -11.
