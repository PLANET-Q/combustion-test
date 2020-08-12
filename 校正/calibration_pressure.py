import numpy as np
from scipy import interpolate

from bokeh.plotting import figure, show, output_file, reset_output
from bokeh.palettes import d3


# 最初に設定するもの----------------------------------------------------------------
# PQ_measure_mainの圧力が下がるときのデータは校正式の算出の時のエラーの原因なので省く
# 最大まで圧力をかけてなだらかに下がっているところまでは入れてもよい
# EDXはしなくてよい
# 自動で除去できるようにしたい
# 測定時に下がっている部分のデータを取らないようにする
filename_pqmain = "PQ_measure_main_36.dat"
filename_edx = "EDX100A_data0182.dat"
calibration_date = "20200114"
# 測定は複数あるので名前を付ける
instrument = "valve"
# 0から数える
# 今後は測定時にどちらを使うかを固定する
pqmain_row = 1
# ----------------------------------------------------------------------------------


# PQ_measure_mainデータの読み込み
folder_filename_pqmain = "input/" + calibration_date + "/" + filename_pqmain
pqmaindata = np.loadtxt(folder_filename_pqmain)
pqmaintime_sec = pqmaindata[:, 0] / 1000000
pqmainpressure = pqmaindata[:, pqmain_row]

# EDXデータの読み込み
folder_filename_edx = "input/" + calibration_date + "/" + filename_edx
edxdata = np.loadtxt(folder_filename_edx)
edxtime = edxdata[:, 0]
edxpressure = edxdata[:, 2]

# データの時間軸をそろえる----------------------------------------------------------
# とんがっているところの先端（最大値）に合わせている
pqmain_pressure_max_time = pqmaintime_sec[pqmainpressure.argmax()]
edx_pressure_max_time = edxtime[edxpressure.argmax()]

if edx_pressure_max_time >= pqmain_pressure_max_time:
    edx_before_len = len(edxtime)
    edxtime = edxtime - (edx_pressure_max_time - pqmain_pressure_max_time)
    edxtime = edxtime[edxtime >= 0]
    edxpressure = edxpressure[(edx_before_len - len(edxtime)):]

if edx_pressure_max_time < pqmain_pressure_max_time:
    pqmain_before_len = len(pqmaintime_sec)
    pqmaintime_sec = pqmaintime_sec - (pqmain_pressure_max_time - edx_pressure_max_time)
    pqmaintime_sec = pqmaintime_sec[pqmaintime_sec >= 0]
    pqmainpressure = pqmainpressure[pqmain_before_len - len(pqmaintime_sec):]

# ----------------------------------------------------------------------------------


# 線形補間
f = interpolate.interp1d(pqmaintime_sec, pqmainpressure, kind = 'linear')

# データの範囲をそろえる
edxtime0_to_pqmaintimesecmin_index = len(edxtime[(edxtime <= pqmaintime_sec.min())])
edxtime = edxtime[(edxtime > pqmaintime_sec.min()) & (edxtime < pqmaintime_sec.max())]
pqmainpressure_linear = f(edxtime)
edxpressure = edxpressure[edxtime0_to_pqmaintimesecmin_index : edxtime0_to_pqmaintimesecmin_index + len(edxtime)]

# 校正係数（一次関数）の算出
coefficient = np.polyfit(pqmainpressure_linear, edxpressure, 1)

pqmainpressure_pma = pqmainpressure_linear * coefficient[0] + coefficient[1]

# 保存
pres_ex_data = np.c_[edxtime, pqmainpressure_pma, edxpressure]
filename_pres_ex = 'output/' + calibration_date + '/pres_ex' + instrument +'.txt'
np.savetxt(filename_pres_ex, pres_ex_data)
fileneme_coefficient = 'output/' + calibration_date + '/coefficient' + instrument + '.txt'
np.savetxt(fileneme_coefficient, coefficient)


# グラフ作成------------------------------------------------------------------------------
reset_output()
filename_result = 'output/' + calibration_date + '/確認' + instrument + '.html'
output_file(filename_result)
C = d3["Category10"][10]

TOOLTIPS_pressure = [("(PQ_measure_main_data, EDX100A_data)", "($x, $y)"), ]
p_pressure = figure(plot_width = 1200, 
                plot_height = 700, 
                x_axis_label = "PQ_measure_main Pressure[MPa]", 
                y_axis_label = "EDX Pressure[MPa]", 
                tooltips = TOOLTIPS_pressure)
p_pressure.line(pqmainpressure_pma, edxpressure, line_color = C[0], 
                line_width = 0.5, legend = "Pressure")

# グリッド
p_pressure.ygrid.minor_grid_line_color = 'lightgray'
p_pressure.ygrid.minor_grid_line_alpha = 0.5
p_pressure.ygrid.minor_grid_line_dash = 'dotted'
p_pressure.xgrid.minor_grid_line_color = 'lightgray'
p_pressure.xgrid.minor_grid_line_dash = 'dotted'
p_pressure.xgrid.minor_grid_line_alpha = 0.5

show(p_pressure)