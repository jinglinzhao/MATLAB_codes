PRO Halpha_Orders,star,no_write=no_write

;Code to compile orders 30 and 31 of UCLES spectra for Jinglin's H-alpha
;measurement pipeline

;Inputs
;star -- string -- star name without the 'HD', e.g. '20807'

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print,' '
print,'Now starting HD '+star
print,' '

;cmdln='mkdir /Users/thedude/idle/universal/jj8/jenns_stuff/LaliotisData/AAT_Work/'+star
;spawn,cmdln

readcol,'AAT_SysRVs.txt',f='A,F',starname,systemicRVs
sysRV=systemicRVs[where(starname eq star)]
;help,sysRV

restore,'/Users/thedude/idle/universal/cfaat/vst'+star+'.dat'
restore,'/Users/thedude/idle/aat/udop/ipcf.dat'
inb=where(ipcf.iodnm eq 'vdiod_' and ipcf.mdchi lt 2.3 and ipcf.cts lt 400000)
qcf=ipcf(inb)

for i=0,n_elements(cf1.obnm)-1 do begin
  spec_exist=findfile('/mir2/iodspec/'+cf1[i].obnm)
  if spec_exist ne '' then begin
    rdsi,spec,cf1[i].obnm 
  endif else begin
    print,'stellar spectrum missing: '+cf1[i].obnm
    print,'skipping this observation'
    goto,skip
  endelse

  cmdln="grep '"+cf1[i].obnm+" "+"' /mir2/bary/ubcvel.ascii"
  spawn,cmdln,output
;  print,output
  targ=strtrim(getwrd(output[0],1),1)
  bc=strtrim(getwrd(output[0],2),1)
  tob=strtrim(getwrd(output[0],3),2)+2440000.
  tob=tob[0]
  tdif=minloc(abs(qcf.jd-tob),/first)
  iodname=qcf(tdif).obnm
  iodstlen=strlen(iodname)
  iodnm=strmid(iodname,0,iodstlen-2)
  exist_b=findfile('/mir2/files/vdiod_'+iodnm+'.b')
  exist_c=findfile('/mir2/files/vdiod_'+iodnm+'.c')

  if exist_c ne '' then begin
    iod_file=exist_c[0] 
  endif else if exist_b ne '' then begin
    iod_file=exist_b[0]
  endif else if (exist_c eq '' && exist_b eq '') then begin
    print,'No iodine found for observation '+cf1[1].obnm
    print,'skipping this file'
    goto,skip
  endif

;  print,'stellar observation: ',cf1[i].obnm
;  print,'closest iodine: ',iod_file

  restore,iod_file
  wceche,vd,wave,/ucles 
  vactoair,wave

  c=299792458 ;in m/s
  bc=float(bc)

  shift= 1 + (sysRV[0]-bc)/c

  halpha_air=6562.808
  halpha_vac=6564.6
  center_line=6563.25

  if keyword_set(plot) then begin
    plot,wave[*,30]/shift,spec[*,30],xr=[6555,6575],tit=star+' '+cf1[i].obnm
    ;oplot,wave[*,30]/shift,spec[*,30],co=151
    oplot,[halpha_air,halpha_air],[0,3e7],co=151
    oplot,[center_line,center_line],[0,3e7],co=90
    ;wait,.1
  endif

  if keyword_set(no_write) then goto,skip
  openw,lun,'/Users/thedude/idle/universal/jj8/jenns_stuff/LaliotisData/AAT_Work/'+star+'/'+cf1[i].obnm+'.txt',/get_lun
  date_str=strtrim(getwrd(output[0],3),2)
  bc_str=strtrim(string(bc),2)
  printf,lun,'# Time: '+date_str+' '
  printf,lun,'# Target: '+'HD'+star+' '
  printf,lun,'# Barycentric correction (m/s): '+bc_str+' '
  printf,lun,'# Spectrum Type: observation'+' '
  printf,lun,'# wavelength flux order'+' '

  dimensions=size(spec)
  n_pixels=dimensions[1]
  
  for j=0,n_pixels-1 do printf,lun,wave[j,30]/shift,'  ',spec[j,30],'  ','30',f='(D15.4,a2,D15.4,a2,a2)'
  for k=0,n_pixels-1 do printf,lun,wave[k,31]/shift,'  ',spec[k,31],'  ','31',f='(D15.4,a2,D15.4,a2,a2)'

  free_lun,lun

  skip:

endfor


END
