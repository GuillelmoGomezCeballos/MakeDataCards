#!/bin/tcsh -f

source ~/EVAL65 5_2_8

#setenv TAG FF_
#setenv TAG SM4_
setenv TAG ""

foreach mass (110 110.5 111 111.5 112 112.5 113 113.5 114 114.5 115 115.5 116 116.5 117 117.5 118 118.5 119 119.5 120 120.5 121 121.5 122 122.5 123 123.5 124 124.1 124.2 124.3 124.4 124.5 124.6 124.7 124.8 124.9 125 125.1 125.2 125.3 125.4 125.5 125.6 125.7 125.8 125.9 126 126.1 126.2 126.3 126.4 126.5 126.6 126.7 126.8 126.9 127 127.5 128 128.5 129 129.5 130 130.5 131 131.5 132 132.5 133 133.5 134 134.5 135 135.5 136 136.5 137 137.5 138 138.5 139 139.5 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 162 164 166 168 170 172 174 176 178 180 182 184 186 188 190 192 194 196 198 200 202 204 206 208 210 212 214 216 218 220 222 224 226 228 230 232 234 236 238 240 242 244 246 248 250 252 254 256 258 260 262 264 266 268 270 272 274 276 278 280 282 284 286 288 290 295 300 305 310 315 320 325 330 335 340 345 350 360 370 380 390 400 420 440 450 460 480 500 520 540 550 560 580 600)

echo $mass

if ( "$TAG" == "FF_" ) then
  sed -ie 's/hwwof_0j.input_7TeV.root/FF_hwwof_0j.input_7TeV.root/' test/$mass/hwwof_0j_shape_7TeV.txt
  #sed -ie 's/hwwsf_0j.input_7TeV.root/FF_hwwsf_0j.input_7TeV.root/' test/$mass/hwwsf_0j_shape_7TeV.txt
  sed -ie 's/hwwof_1j.input_7TeV.root/FF_hwwof_1j.input_7TeV.root/' test/$mass/hwwof_1j_shape_7TeV.txt
  #sed -ie 's/hwwsf_1j.input_7TeV.root/FF_hwwsf_1j.input_7TeV.root/' test/$mass/hwwsf_1j_shape_7TeV.txt
endif
if ( "$TAG" == "SM4_" ) then
  sed -ie 's/hwwof_0j.input_7TeV.root/SM4_hwwof_0j.input_7TeV.root/' test/$mass/hwwof_0j_shape_7TeV.txt
  #sed -ie 's/hwwsf_0j.input_7TeV.root/SM4_hwwsf_0j.input_7TeV.root/' test/$mass/hwwsf_0j_shape_7TeV.txt
  sed -ie 's/hwwof_1j.input_7TeV.root/SM4_hwwof_1j.input_7TeV.root/' test/$mass/hwwof_1j_shape_7TeV.txt
  #sed -ie 's/hwwsf_1j.input_7TeV.root/SM4_hwwsf_1j.input_7TeV.root/' test/$mass/hwwsf_1j_shape_7TeV.txt
endif

mkdir -p /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass

cp test/$mass/hwwof_0j_shape_7TeV.txt   /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwof_0j_shape_7TeV.txt 
#cp test/$mass/hwwsf_0j_shape_7TeV.txt   /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwsf_0j_shape_7TeV.txt 
cp test/$mass/hwwof_1j_shape_7TeV.txt   /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwof_1j_shape_7TeV.txt 
#cp test/$mass/hwwsf_1j_shape_7TeV.txt   /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwsf_1j_shape_7TeV.txt 
#cp test/$mass/hwwof_2j_shape_7TeV.txt   /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwof_2j_shape_7TeV.txt 
#cp test/$mass/hwwsf_2j_shape_7TeV.txt   /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwsf_2j_shape_7TeV.txt 
cp test/$mass/hwwof_0j.input_7TeV.root  /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwof_0j.input_7TeV.root
#cp test/$mass/hwwsf_0j.input_7TeV.root  /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwsf_0j.input_7TeV.root
cp test/$mass/hwwof_1j.input_7TeV.root  /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwof_1j.input_7TeV.root
#cp test/$mass/hwwsf_1j.input_7TeV.root  /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwsf_1j.input_7TeV.root
#cp test/$mass/hwwof_2j.input_7TeV.root  /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwof_2j.input_7TeV.root
#cp test/$mass/hwwsf_2j.input_7TeV.root  /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwsf_2j.input_7TeV.root
cp test/$mass/hwwof_0j_cut_7TeV.txt	/build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwof_0j_cut_7TeV.txt
cp test/$mass/hwwsf_0j_cut_7TeV.txt	/build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwsf_0j_cut_7TeV.txt
cp test/$mass/hwwof_1j_cut_7TeV.txt	/build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwof_1j_cut_7TeV.txt
cp test/$mass/hwwsf_1j_cut_7TeV.txt	/build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwsf_1j_cut_7TeV.txt
cp test/$mass/hwwof_2j_cut_7TeV.txt     /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwof_2j_cut_7TeV.txt   
cp test/$mass/hwwsf_2j_cut_7TeV.txt     /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwsf_2j_cut_7TeV.txt   

#text2workspace.py -b /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwof_0j_shape_7TeV.txt 
#text2workspace.py -b /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwsf_0j_shape_7TeV.txt 
#rm -f /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwof_0j_shape_7TeV.root
#rm -f /build/ceballos/cmshcg/trunk/summer2013/couplings/hww2l2v/$mass/${TAG}hwwsf_0j_shape_7TeV.root

end
