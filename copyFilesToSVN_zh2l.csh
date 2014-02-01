#!/bin/tcsh -f

source ~/EVAL65 5_2_8

#setenv TAG FF_
#setenv TAG SM4_
setenv TAG ""

setenv MAINDIR /build/ceballos/cmshcg/trunk/summer2013/couplings/zllhinv

#couplings
#foreach mass (110 110.5 111 111.5 112 112.5 113 113.5 114 114.5 115 115.5 116 116.5 117 117.5 118 118.5 119 119.5 120 120.5 121 121.5 122 122.5 123 123.5 124 124.1 124.2 124.3 124.4 124.5 124.6 124.7 124.8 124.9 125 125.1 125.2 125.3 125.4 125.5 125.6 125.7 125.8 125.9 126 126.1 126.2 126.3 126.4 126.5 126.6 126.7 126.8 126.9 127 127.5 128 128.5 129 129.5 130 130.5 131 131.5 132 132.5 133 133.5 134 134.5 135 135.5 136 136.5 137 137.5 138 138.5 139 139.5 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 162 164 166 168 170 172 174 176 178 180 182 184 186 188 190 192 194 196 198 200 300)
foreach mass (115 135 145 200 300 124.5 124.6 124.7 124.8 124.9 125 125.1 125.2 125.3 125.4 125.5 125.6 125.7 125.8 125.9 126 126.1 126.2 126.3 126.4 126.5)

echo $mass

mkdir -p ${MAINDIR}/$mass

cp test/$mass/zmmhinv_nj0_shape_7TeV.txt    ${MAINDIR}/$mass/${TAG}zmmhinv_nj0_shape_7TeV.txt
cp test/$mass/zeehinv_nj0_shape_7TeV.txt    ${MAINDIR}/$mass/${TAG}zeehinv_nj0_shape_7TeV.txt
cp test/$mass/zmmhinv_nj1_shape_7TeV.txt    ${MAINDIR}/$mass/${TAG}zmmhinv_nj1_shape_7TeV.txt
cp test/$mass/zeehinv_nj1_shape_7TeV.txt    ${MAINDIR}/$mass/${TAG}zeehinv_nj1_shape_7TeV.txt

cp test/$mass/zllhinvmm0j_input_7TeV.root    ${MAINDIR}/$mass/${TAG}zllhinvmm0j_input_7TeV.root
cp test/$mass/zllhinvee0j_input_7TeV.root    ${MAINDIR}/$mass/${TAG}zllhinvee0j_input_7TeV.root
cp test/$mass/zllhinvmm1j_input_7TeV.root    ${MAINDIR}/$mass/${TAG}zllhinvmm1j_input_7TeV.root
cp test/$mass/zllhinvee1j_input_7TeV.root    ${MAINDIR}/$mass/${TAG}zllhinvee1j_input_7TeV.root

end

foreach mass (105 175)

echo $mass

mkdir -p ${MAINDIR}/$mass

cp /data/smurf/ceballos/inputLimits/ana_ZHinv/$mass/zmmhinv_nj0_shape_7TeV.txt	${MAINDIR}/$mass/${TAG}zmmhinv_nj0_shape_7TeV.txt
cp /data/smurf/ceballos/inputLimits/ana_ZHinv/$mass/zeehinv_nj0_shape_7TeV.txt	${MAINDIR}/$mass/${TAG}zeehinv_nj0_shape_7TeV.txt
cp /data/smurf/ceballos/inputLimits/ana_ZHinv/$mass/zmmhinv_nj1_shape_7TeV.txt	${MAINDIR}/$mass/${TAG}zmmhinv_nj1_shape_7TeV.txt
cp /data/smurf/ceballos/inputLimits/ana_ZHinv/$mass/zeehinv_nj1_shape_7TeV.txt	${MAINDIR}/$mass/${TAG}zeehinv_nj1_shape_7TeV.txt

cp /data/smurf/ceballos/inputLimits/ana_ZHinv/$mass/zllhinvmm0j_input_7TeV.root	 ${MAINDIR}/$mass/${TAG}zllhinvmm0j_input_7TeV.root
cp /data/smurf/ceballos/inputLimits/ana_ZHinv/$mass/zllhinvee0j_input_7TeV.root	 ${MAINDIR}/$mass/${TAG}zllhinvee0j_input_7TeV.root
cp /data/smurf/ceballos/inputLimits/ana_ZHinv/$mass/zllhinvmm1j_input_7TeV.root	 ${MAINDIR}/$mass/${TAG}zllhinvmm1j_input_7TeV.root
cp /data/smurf/ceballos/inputLimits/ana_ZHinv/$mass/zllhinvee1j_input_7TeV.root	 ${MAINDIR}/$mass/${TAG}zllhinvee1j_input_7TeV.root

end
