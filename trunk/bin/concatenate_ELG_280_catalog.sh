#!/bin/bash/

ls ../decalsSelections/bricks/*ELG_280*.fits > flAll
sed -n '1,10000p' flAll > flAll_1
sed -n '10001,20000p' flAll > flAll_2

java -jar stilts.jar tcat ifmt=fits in=@flAll_1 out=../decalsSelections/catalog_ELG_280_1.fits
java -jar stilts.jar tcat ifmt=fits in=@flAll_2 out=../decalsSelections/catalog_ELG_280_2.fits

ls ../decalsSelections/*ELG_280_?.fits > flAll_fits
java -jar stilts.jar tcat ifmt=fits in=@flAll_fits out=../decalsSelections/catalog_ELG_280.fits

rm flAll*
rm ../decalsSelections/catalog_ELG_2?0_?.fits

