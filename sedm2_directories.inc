;----------- directories info------
  ;dir the base directory that store the filters and SSP models
  dir = '/home/yz69/SEDMORPH/'

  ; SDSS filters, ask V. Wild for more info
  dir_filters = dir+'FILTERS/'

  ; SSP models, here we use BC03 models with Chabrier IMF
  ; check http://www.bruzual.org/
  ; dir_models = dir+'Miles_Atlas/Chabrier_IMF/'
  ; We interpolated a few metallicity bins, check Yirui Zheng+2020 for more details 
  dir_models = dir+'Miles_Atlas/Chabrier_IMF_interpolated/'

  ; Where you store the PSF weight files. The files can be created with
  ; .....ipynb in the subdirectory "some_dependencies"
  ; check Yirui Zheng+2020 and his upcoming thesis for more info
  dir_PSF_weight = dir+'public_SEDMORPH_SIMN2/some_dependencies/PSF_mass_weight/'

  ; dir_IN, where you store the simulation data, i.e. the snapshot files
  dir_IN = dir+'SimnsGadget3/data/major/'+fileseq+'/'
  ;dir_IN = dir+'SimnsGadget3/data/minor/'+fileseq+'/'
  ;dir_IN = dir+'SimnsGadget3/data/isolated/'+fileseq+'/'

  ; dir_OUT where you store the output image, spectra, pcs and so on
  ; dir_OUT = dir+'SimnsGadget3_output/'+fileseq+'/'
  dir_OUT = '/share/data/yz69/SEDMORPH/'+'SimnsGadget3_output/'+fileseq+'/'

  dir_code = '/home/yz69/SEDMORPH/public_SEDMORPH_SIMN2' ;; needed for galaxy image which is stored here

  ; where you store the pca data
  ; Check http://www-star.st-and.ac.uk/~vw8/downloads/DR7PCA.html
  ; cite  Wild et al. 2007, MNRAS, 381, 543 
  dir_pca_data = dir+'pca_data/'
