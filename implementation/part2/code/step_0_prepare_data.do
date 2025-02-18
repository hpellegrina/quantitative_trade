clear all

* make microregion mask
use ".././input/data/geo_dist_mat.dta", clear
keep micro_code_i
rename micro_code_i micro_code
duplicates drop micro_code, force
save "../output/microregion_mask.dta", replace

sort micro_code 
export delimited micro_code using ".././output/micro_code.csv", novar replace

* bring population
import excel ".././input/data/pop.xls", sheet("Séries") firstrow clear
keep Codigo D
rename (Codigo D) (micro_code pop)
destring micro_code, replace
merge 1:1 micro_code using "../output/microregion_mask.dta", nogen keep(3)
sort micro_code 
export delimited pop using ".././output/pop.csv", novar replace

* bring gdp
import excel ".././input/data/gdp.xls", sheet("Séries") firstrow clear
keep Codigo D
rename (Codigo D) (micro_code gdp)
destring micro_code, replace
merge 1:1 micro_code using "../output/microregion_mask.dta", nogen keep(3)
sort micro_code 
export delimited gdp using ".././output/gdp.csv", novar replace

* export distances
use ".././input/data/geo_dist_mat.dta", clear
keep micro_code_i micro_code_j dist
sort micro_code_i micro_code_j
reshape wide dist, i(micro_code_i) j(micro_code_j)
drop micro_code_i
export delimited dist* using ".././output/mat_dist.csv", novar replace

* erase temporary files
erase "../output/microregion_mask.dta"
