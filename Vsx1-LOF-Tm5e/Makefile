.PHONY: clean
INPUTS := $(wildcard ./batch1/mCherryi_*_D_Ey_Kn.tif ./batch1/VKD_*_D_Ey_Kn.tif ./batch2/mCherryi_*_D_Ey_Kn.tif ./batch2/VKD_*_D_Ey_Kn.tif)

.binarize.done: $(INPUTS)
	bash ./get_dbp.sh
	touch .binarize.done

binarize: .binarize.done

clean:
	find ./ -name "*eykn*" -delete
	rm .binarize.done
