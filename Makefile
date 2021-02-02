all:
	cd util/misc/data && ./build_rg_obj.R


clean:
	cd util/misc/data && rm -f ./ranger.rg_obj.rds
