;;*** CSV_STR.PRO
;;
;; AIM: turn a csv string array into a structure array
;;
;; INPUT: 
;;        tags: tag names
;;        csv_result: string array
;; OPTIONAL INPUT: 
;;        format: format of string (necessary if there might be a gap
;;        in the first entry, or a confusion between long and int),
;;        same as readcol string of 'A,F,I,L,D'
;;
;; OUTPUT:
;;        result: structure array
;; 
;;------------------------------------------------------------------

FUNCTION csv_str, tags,csv_result,format=format,silent=silent



;sample data:

;	# Result of SQL query :
;	#select * from agebin
;	runId,binId,minAge,maxAge,cenAge,minRedshift,maxRedshift,cenRedshift
;	13,0,0.00965843,0.02,0.0138985,7.016178E-4,0.0014536649,0.0010098685
;	13,1,0.02,0.0414144,0.02878,0.0014536649,0.0030134656,0.0020927698


fieldnames = STRSPLIT(tags, ',', /extract)
n_fields=n_elements(fieldnames)
n_records=n_elements(csv_result)

;get rid of laeding and trailing blanks
for i=0,n_fields-1 do fieldnames[i]=strtrim(fieldnames[i],2)

;create structure
q_fn="'"+fieldnames+"'"

;; vw: loop added incase of null fields in first row. This needs improving!
if n_elements(format) eq 0 then begin
    i=0L
    while (1) do begin
        sample_data = strsplit(csv_result[i], ',', /extract)
        if n_elements(sample_data) eq n_fields then break else i=i+1L
        if i eq n_records -1 then begin
            print,'no complete data row: please use format keyword'
            return,0
        endif
    endwhile
    
;check for spaces or letters in data, which implies a string
;first remove leading and trailing white spaces
    for i=0,n_fields-1 do sample_data[i]=strtrim(sample_data[i],2)
;put quotations around string
    for i=0, n_fields-1 do begin
        if (strpos(sample_data[i], ' ') ne -1 OR $ ;does it have a space in?
            STREGEX(sample_data[i], '[a,b,c,d,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z]', /FOLD_CASE) ne -1) AND $ ;does it have a letter in? (e could be nb, not ideal)
          strmid(sample_data[i], 0, 1) ne "'" then $ ;is it already within quotes?
          sample_data[i]="'"+sample_data[i]+"'" 
    endfor
endif else begin
    format_str = strsplit(format,',',/extract)
    sample_data = strarr(n_fields)
    if n_elements(format_str) ne n_fields then begin
        print,'!!CSV_STR: size of format not equal to size of data, returning'
        return,-1
    endif
    for i=0,n_fields-1 do begin
        if strmatch(format_str[i],'*A*') then sample_data[i] = "'a'"
        if strmatch(format_str[i],'*F*') then sample_data[i] = '0.0'
        if strmatch(format_str[i],'*D*') then sample_data[i] = '0D'
        if strmatch(format_str[i],'*I*') then sample_data[i] = '0'
        if strmatch(format_str[i],'*L*') then sample_data[i] = '0L'
        if strmatch(format_str[i],'*V*') then sample_data[i] = 'long64(0)'
    endfor
endelse




big_box_str=transpose([[q_fn], [sample_data]])
big_line_str=strjoin(strjoin(big_box_str, ', '), ', ')

comm='res=create_struct('+big_line_str+')'
ro=execute(comm)

if not(keyword_set(silent)) then print, 'Parsing a search returning ' + string(n_records) + ' lines'

result=replicate(res, n_records)

for i=0l, n_records-1 do begin
    vals=strsplit(csv_result[i], ',', /extract,/preserve_null)
    ind = where(vals eq 'null')
    if ind[0] ne -1 then vals[ind]='-999'

    for j=0,n_fields-1 do result[i].(j)=strtrim(vals[j],2)
end

return, result
	
END
