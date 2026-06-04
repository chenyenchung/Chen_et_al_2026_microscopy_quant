.PHONY: clean
INPUTS := $(wildcard ./batch1/LacZ_*_D_Ey_Kn.tif ./batch1/Vsx1_*_D_Ey_Kn.tif ./batch2/LacZ_*_D_Ey_Kn.tif ./batch2/Vsx1_*_D_Ey_Kn.tif)

.binarize.done: $(INPUTS)
	bash ./get_dbp.sh
	touch .binarize.done

binarize: .binarize.done

clean:
	find ./ -name "*eykn*" -delete
	find ./ -name "*cp_mask*" -delete
	rm .binarize.done
