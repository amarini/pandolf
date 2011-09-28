# Simple counting experiment, with one signal and one background process
#imax 1  number of channels
#jmax 1  number of backgrounds 
#kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------
shapes ggH CMS_hwwlvqq_mu hwwlvqq_mu.input.root  w:signal
shapes VBF CMS_hwwlvqq_mu hwwlvqq_mu.input.root  w:signal 
shapes background CMS_hwwlvqq_mu hwwlvqq_mu.input.root w:background
###shapes data_obs CMS_hwwlvqq_mu hwwlvqq_mu.input.root w:data_obs
shapes data_obs CMS_hwwlvqq_mu hwwlvqq_mu.input.root w:dataset_obs
------------
bin         CMS_hwwlvqq_mu
observation <dummyobs>
------------
bin                CMS_hwwlvqq_mu	  CMS_hwwlvqq_mu	          CMS_hwwlvqq_mu
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
CMS_trigger_mu	lnN	1.02	1.02	1.0	
CMS_eff_mu	lnN	1.008	1.008	1.0
CMS_scale_mu	lnN	1.01	1.01	1.0
#CMS_reco_mu	lnN	1.015	1.015	1.0
CMS_scale_j	<dummyJES>
CMS_eff_b 	<dummybeff>
CMS_hwwlvqq_pu	lnN		1.02	        1.02		1.0			      	
CMS_hwwlvqq_qgsep   lnN	1.046            1.046          1.0
#CMS_hwwlvqq_sigp0
#...                                                  <developing...>
#CMS_hwwlvqq_sigp5
CMS_hwwlvqq_bkgmup0   <dummybnorm>
#CMS_hwwlvqq_bkgp1        param  222.72  ---
CMS_hwwlvqq_bkgp2         param -0.116428   0.0453989
#CMS_hwwlvqq_bkgp3      param   13.3992   ----
CMS_hwwlvqq_bkgp4         param 54.7039  6.792

