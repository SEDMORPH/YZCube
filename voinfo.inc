
;; Information for creating eigenspectra
DIR_VOpca = dir_pca_data+'VO/PCARUN/'
info_pca = {runno:0, time:'',data:'',description:'',meansub:0,varweight:0,EMpca:0,EMpca_niter:0,gappy:0,gappy_niter:0,vdisp:0.0,norm:'',minwave:0.,maxwave:0.,extra:''}
VOinfofile = DIR_VOpca+'pcavo_info.csv'

;; Information for projecting eigenspectra onto different catalogues
DIR_VOproj = '/scratch/vwild/VO/PROJECT/'
info_proj = {projno:0,especno:0,time:'',catalogue:'',selection:'', gappy:0, dust:0}
PROJinfofile = DIR_VOproj+'pcavo_proj.csv'
