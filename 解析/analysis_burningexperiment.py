import numpy as np
from scipy import fftpack
import datetime

from bokeh.plotting import figure, save, output_file, reset_output
from bokeh.palettes import d3
from bokeh.models.annotations import Label
from bokeh.layouts import gridplot



#最初に設定する値---------------------------------------------------------------------------------------------------------------------------
filename_edx = "EDX100A_data0188.dat"
filename_pqmain = "PQ_measure_main_51.dat"

experiment_date = "20200211"
enginename = "ntk360_OC500"

#温度圧力の順番、名前も必要に応じて変える
#PQ_measure_mainの2行目　燃焼室圧
#                 3     タンク圧
#     　　　　　　 4     タンク温度　　　
#　　   　　　　　 5        |
#　　　 　  　　　 6        |
#                 7　　 下部燃焼室温度
#今後は圧力と温度の順番を固定する
temperature1_name = "lower pipe"
temperature2_name = "upper engine"
temperature3_name = "tank"
tempereture4_name = "upper pipe"

#校正係数傾き
inclination_calibration_coefficient_tank = 0.000654202
inclination_calibration_coefficient_firing = 0.000654202
#校正係数y切片
intercept_calibration_coefficient_tank = -1.733407
intercept_calibration_coefficient_firing = -1.733407

#ノイズを除去したグラフが変だったら、filename_freqs_power_thrustのグラフを確認して、カットオフ周波数を変えてみると直る
#0.0000002で大体のデータはOKっぽい
cutoff_frequency = 0.0000002

#設定する値はここまで-----------------------------------------------------------------------------------------------------------------------



#EDXデータの読み込み------------------------------------------------------------------------------------------------------------------------
folder_filename_edx = "input/" + experiment_date + "/" + filename_edx
edxdata = np.loadtxt(folder_filename_edx, usecols = (0, 1))

#0合わせ
edxdata_offset_select = edxdata[100:200, 1]
edx_offset_ave = np.mean(edxdata_offset_select)
thrust = edxdata[:, 1] - edx_offset_ave

edxtime = edxdata[:, 0]
#-------------------------------------------------------------------------------------------------------------------------------------------


#ノイズ除去---------------------------------------------------------------------------------------------------------------------------------
#カットオフ周波数を求める
#今はカットオフ周波数のグラフから自動でカットオフ周波数を求める方法がわからない
#0.0000002で大体のデータはOKっぽい
sample_freq_thrust =fftpack.fftfreq(thrust.size,d=10000)
sig_fft_thrust=fftpack.fft(thrust)
pidxs_thrust = np.where(sample_freq_thrust > 0)
freqs_thrust, power_thrust = sample_freq_thrust[pidxs_thrust], np.abs(sig_fft_thrust)[pidxs_thrust]
freqs_power_thrust = np.c_[freqs_thrust, power_thrust]

#保存する
now = datetime.datetime.now()
filename_freqs_power_thrust = 'output/' + experiment_date + '/freqs_power_thrust' + now.strftime('%Y%m%d_%H%M%S') + '.txt'
np.savetxt(filename_freqs_power_thrust, freqs_power_thrust)

#FFT計算
sig_fft_thrust[np.abs(sample_freq_thrust) > cutoff_frequency] = 0
fft_thrust =np.real( fftpack.ifft(sig_fft_thrust))

#FFT保存する
now = datetime.datetime.now()
filename_fft_thrust = 'output/' + experiment_date + '/fft_thrust' + now.strftime('%Y%m%d_%H%M%S') + '.txt'
np.savetxt(filename_fft_thrust, np.c_[edxtime, fft_thrust])
#-------------------------------------------------------------------------------------------------------------------------------------------


#燃焼時間計算-------------------------------------------------------------------------------------------------------------------------------
#作動開始、作動終了のindexの算出-------------------------------------------------------------------
#最大推力の10%となる推力を計算している
thrust_max = fft_thrust.max()
thrust_max_index = fft_thrust.argmax()
thrust_10per = thrust_max * 0.1
thrust_aftermax = fft_thrust[thrust_max_index:]
thrust_beforemax = fft_thrust[:thrust_max_index]
thrust_aftermax_large10per = thrust_aftermax[thrust_aftermax > thrust_10per]
action_last = len(thrust_aftermax_large10per) + thrust_max_index
thrust_beforemax_small10per = thrust_beforemax[thrust_beforemax < thrust_10per]
action_first = len(thrust_beforemax_small10per)

#燃焼終了のindexの算出-----------------------------------------------------------------------------
#ノイズ除去推力の読み込み
fft_data = np.loadtxt(filename_fft_thrust)
fft_data_action = fft_data[action_first:action_last, :]

#燃焼終了を求めるために各点で傾きとy切片を求めている
differential = (fft_data_action[1:, 1] - fft_data_action[:-1, 1]) / (fft_data_action[1:, 0] - fft_data_action[:-1, 0])
y_intercept = differential * (- fft_data_action[1:,0]) + fft_data_action[1:, 1]

#傾きとy切片の保存
now = datetime.datetime.now()
filename_tangent_line = 'output/' + experiment_date + '/tangent_line' + now.strftime('%Y%m%d_%H%M%S') + '.txt'
np.savetxt(filename_tangent_line, np.c_[differential, y_intercept])

#y切片が最大の点からindexを小さくしていって、最初にy切片が0となる点のindexを求めて燃焼終了のindexとする
#TEKだとうまくいかない
y_intercept_max_index = y_intercept.argmax()
y_intercept_beforemax = y_intercept[:y_intercept_max_index]
count_y_intercept_beforemax_to_0 = 0
for c in reversed(y_intercept_beforemax):
    if c < 0:
        break
    count_y_intercept_beforemax_to_0 = count_y_intercept_beforemax_to_0 + 1
burning_last = y_intercept_max_index - count_y_intercept_beforemax_to_0 + action_first

#燃焼時間と作動時間の算出--------------------------------------------------------------------------
burningtime = edxtime[burning_last] - edxtime[action_first]
actiontime = edxtime[action_last] - edxtime[action_first]

#燃焼時間計算終了---------------------------------------------------------------------------------------------------------------------------


#平均推力計算
actiontime_thrust_ave = np.mean(fft_thrust[action_first:action_last])
burningtime_thrust_ave = np.mean(fft_thrust[action_first:burning_last])

#Total Impulse計算
impulse = (edxtime[1:] - edxtime[:-1]) * (fft_thrust[1:] + fft_thrust[:-1]) / 2
total_impulse = np.sum(impulse[action_first:action_last])


#自作測定器のデータ読み込み-----------------------------------------------------------------------------------------------------------------
folder_filename_pqmain = "input/" + experiment_date + "/" + filename_pqmain
pqmaindata = np.loadtxt(folder_filename_pqmain)
pqmaintime_sec = pqmaindata[:, 0] / 1000000
firingpressure = pqmaindata[:, 1] * inclination_calibration_coefficient_firing + intercept_calibration_coefficient_firing
tankpressure = pqmaindata[:, 2] * inclination_calibration_coefficient_tank + intercept_calibration_coefficient_tank
temperature1 = pqmaindata[:, 3]
temperature2 = pqmaindata[:, 4]
temperature3 = pqmaindata[:, 5]
temperature4 = pqmaindata[:, 6]
#-------------------------------------------------------------------------------------------------------------------------------------------



#グラフ作成---------------------------------------------------------------------------------------------------------------------------------
reset_output()
filename_result = 'output/' + experiment_date + '/' + enginename + '_燃焼試験解析.html'
output_file(filename_result)
C = d3["Category10"][10]

#グリッドを統一して整理したほうがいい
#推力----------------------------------------------------------------------------------------------
TOOLTIPS_thrust = [("(time[s], Thrust[N])", "($x, $y)"), ]
p_thrust = figure(title = "thrust", 
                  plot_width = 1200, 
                  plot_height = 700, 
                  x_axis_label = "Time[s]", 
                  y_axis_label = "Thrust[N]", 
                  tooltips = TOOLTIPS_thrust)
p_thrust.line(edxtime, thrust, line_color = C[0], 
              line_width = 0.5, legend = "Thrust")
p_thrust.line(edxtime, fft_thrust, 
              line_color = C[1], legend = "smoothed_data")
p_thrust.line([0, edxtime.max()], thrust_max, 
              line_color = C[2], legend = "max")
p_thrust.line([0, edxtime.max()], thrust_10per, 
              line_color = C[3], legend = "10percent")
p_thrust.line([0, edxtime.max()], actiontime_thrust_ave, 
              line_color = C[4], legend = "actiontime_average")
p_thrust.line([0, edxtime.max()], burningtime_thrust_ave, 
              line_color = C[5], legend = "burningtime_average")
p_thrust.line(edxtime[action_first], [0, thrust.max()], 
              line_color = C[6], legend = "burningstart")
p_thrust.line(edxtime[burning_last], [0, thrust.max()], 
              line_color = C[7], legend = "burningfinish")
p_thrust.line(edxtime[action_last], [0, thrust.max()], 
              line_color = C[8], legend = "actionfinish")
p_thrust.legend.click_policy = "hide"

#ラベル
label_totalimpulse=Label(x=0, y=200, text="total impulse = " + str(total_impulse) + "[Ns]")
p_thrust.add_layout(label_totalimpulse)
label_actiontime=Label(x=0, y=250, text="action time = " + str(actiontime) + "[s]")
p_thrust.add_layout(label_actiontime)
label_burningtime=Label(x=0, y=300, text="burning time = " + str(burningtime) + "[s]")
p_thrust.add_layout(label_burningtime)

#グリッド
p_thrust.ygrid.minor_grid_line_color = 'lightgray'
p_thrust.ygrid.minor_grid_line_alpha = 0.5 #透明度
p_thrust.ygrid.minor_grid_line_dash = 'dotted'
p_thrust.xgrid.minor_grid_line_color = 'lightgray'
p_thrust.xgrid.minor_grid_line_color = 'lightgray'
p_thrust.xgrid.minor_grid_line_alpha = 0.5


#タンク圧------------------------------------------------------------------------------------------
TOOLTIPS_tank = [("(time[s], Tank Pressure[MPa])", "($x, $y)"), ]
p_tank = figure(title = "Tank Pressure", 
                plot_width = 1200, 
                plot_height = 700, 
                x_axis_label = "Time[s]", 
                y_axis_label = "Tank Pressure[MPa]", 
                x_range = p_thrust.x_range, 
                tooltips = TOOLTIPS_tank)
p_tank.line(pqmaintime_sec, tankpressure, line_color = C[0], 
            line_width = 0.5, legend = "Tank Pressure")

#グリッド
p_tank.ygrid.minor_grid_line_color = 'lightgray'
p_tank.ygrid.minor_grid_line_alpha = 0.5
p_tank.ygrid.minor_grid_line_dash = 'dotted'
p_tank.xgrid.minor_grid_line_color = 'lightgray'
p_tank.xgrid.minor_grid_line_dash = 'dotted'
p_tank.xgrid.minor_grid_line_alpha = 0.5


#燃焼室圧------------------------------------------------------------------------------------------
TOOLTIPS_firing = [("(time[s], Firing Pressure[MPa])", "($x, $y)"), ]
p_firing = figure(title = "Firing Pressure", 
                  plot_width = 1200, 
                  plot_height = 700, 
                  x_axis_label = "Time[s]", 
                  y_axis_label = "Firing Pressure[MPa]", 
                  x_range = p_thrust.x_range, 
                  tooltips = TOOLTIPS_firing)
p_firing.line(pqmaintime_sec, firingpressure, line_color = C[1], 
              line_width = 0.5, legend = "Firing Pressure")
#グリッド
p_firing.ygrid.minor_grid_line_color = 'lightgray'
p_firing.ygrid.minor_grid_line_alpha = 0.5
p_firing.ygrid.minor_grid_line_dash = 'dotted'
p_firing.xgrid.minor_grid_line_color = 'lightgray'
p_firing.xgrid.minor_grid_line_dash = 'dotted'
p_firing.xgrid.minor_grid_line_alpha = 0.5

TOOLTIPS_temperature = [("(time[s], Tank Temperature[℃])", "($x, $y)"), ]


#温度1---------------------------------------------------------------------------------------------
p_temperature1 = figure(title = temperature1_name, 
                        plot_width = 700, 
                        plot_height = 500, 
                        x_axis_label = "Time[s]", 
                        y_axis_label = "Temperature[℃]", 
                        x_range = p_thrust.x_range, 
                        tooltips = TOOLTIPS_temperature)
p_temperature1.line(pqmaintime_sec, temperature1, line_color = C[2], 
                    line_width = 0.5, legend = temperature1_name)
#グリッド
p_temperature1.ygrid.minor_grid_line_color = 'lightgray'
p_temperature1.ygrid.minor_grid_line_alpha = 0.5
p_temperature1.ygrid.minor_grid_line_dash = 'dotted'
p_temperature1.xgrid.minor_grid_line_color = 'lightgray'
p_temperature1.xgrid.minor_grid_line_dash = 'dotted'
p_temperature1.xgrid.minor_grid_line_alpha = 0.5


#温度2---------------------------------------------------------------------------------------------
p_temperature2 = figure(title = temperature2_name, 
                        plot_width = 700, 
                        plot_height = 500, 
                        x_axis_label = "Time[s]", 
                        y_axis_label = "Temperature[℃]", 
                        x_range = p_thrust.x_range, 
                        tooltips = TOOLTIPS_temperature)
p_temperature2.line(pqmaintime_sec, temperature2, line_color = C[3], 
                    line_width = 0.5, legend = temperature2_name)
#グリッド
p_temperature2.ygrid.minor_grid_line_color = 'lightgray'
p_temperature2.ygrid.minor_grid_line_alpha = 0.5
p_temperature2.ygrid.minor_grid_line_dash = 'dotted'
p_temperature2.xgrid.minor_grid_line_color = 'lightgray'
p_temperature2.xgrid.minor_grid_line_dash = 'dotted'
p_temperature2.xgrid.minor_grid_line_alpha = 0.5


#温度3---------------------------------------------------------------------------------------------
p_temperature3 = figure(title = temperature3_name, 
                        plot_width = 700, 
                        plot_height = 500, 
                        x_axis_label = "Time[s]", 
                        y_axis_label = "Temperature[℃]", 
                        x_range = p_thrust.x_range, 
                        tooltips = TOOLTIPS_temperature)
p_temperature3.line(pqmaintime_sec, temperature3, line_color = C[4], 
                    line_width = 0.5, legend = temperature3_name)
#グリッド
p_temperature3.ygrid.minor_grid_line_color = 'lightgray'
p_temperature3.ygrid.minor_grid_line_alpha = 0.5
p_temperature3.ygrid.minor_grid_line_dash = 'dotted'
p_temperature3.xgrid.minor_grid_line_color = 'lightgray'
p_temperature3.xgrid.minor_grid_line_dash = 'dotted'
p_temperature3.xgrid.minor_grid_line_alpha = 0.5


#温度4---------------------------------------------------------------------------------------------
p_temperature4 = figure(title = tempereture4_name, 
                        plot_width = 700, 
                        plot_height = 500, 
                        x_axis_label = "Time[s]", 
                        y_axis_label = "Temperature[℃]", 
                        x_range = p_thrust.x_range, 
                        tooltips = TOOLTIPS_temperature)
p_temperature4.line(pqmaintime_sec, temperature4, line_color = C[5], 
                    line_width = 0.5, legend = tempereture4_name)
#グリッド
p_temperature4.ygrid.minor_grid_line_color = 'lightgray'
p_temperature4.ygrid.minor_grid_line_alpha = 0.5
p_temperature4.ygrid.minor_grid_line_dash = 'dotted'
p_temperature4.xgrid.minor_grid_line_color = 'lightgray'
p_temperature4.xgrid.minor_grid_line_dash = 'dotted'
p_temperature4.xgrid.minor_grid_line_alpha = 0.5
#--------------------------------------------------------------------------------------------------

#配置
plot = gridplot([[p_thrust], [p_tank], [p_firing], 
                [p_temperature1, p_temperature2], [p_temperature3, p_temperature4]])
save(plot)
