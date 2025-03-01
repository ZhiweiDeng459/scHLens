import seaborn as sns
from matplotlib.colors import to_hex
import random
import copy

## color_list
color_list = []
color_list.extend(list(sns.color_palette("colorblind").as_hex()))
color_list.extend(list(sns.color_palette("pastel").as_hex()))
color_list.extend(list(sns.color_palette("husl").as_hex()))
color_list.extend(list(sns.color_palette("Set2").as_hex()))
# random_color = [to_hex(sns.color_palette("Spectral", as_cmap=True)(random.random())) for i in range(0,300)]
append_color = ['#9e0142', '#a20643', '#a90d45', '#ad1246', '#b41947', '#b81e48', '#be254a', '#c32a4b', '#c9314c', '#d0384e', '#d43d4f', '#d8434e', '#da464d', '#de4c4b', '#e1504b', '#e45549', '#e75948', '#ea5e47', '#ee6445', '#f06744', '#f46d43', '#f57245', '#f67a49', '#f67f4b', '#f8864f', '#f98e52', '#f99355', '#fa9b58', '#fba05b', '#fca85e', '#fdad60', '#fdb365', '#fdb768', '#fdbd6d', '#fdc372', '#fdc776', '#fecc7b', '#fed07e', '#fed683', '#feda86', '#fee08b', '#fee28f', '#fee695', '#feea9b', '#feec9f', '#fff0a6', '#fff2aa', '#fff6b0', '#fff8b4', '#fffcba', '#ffffbe', '#fdfebb', '#fafdb7', '#f8fcb4', '#f5fbaf', '#f3faac', '#f0f9a7', '#eef8a4', '#ebf7a0', '#e8f69b', '#e6f598', '#dff299', '#daf09a', '#d3ed9c', '#cfec9d', '#c8e99e', '#c3e79f', '#bce4a0', '#b5e1a2', '#b1dfa3', '#aadca4', '#a4daa4', '#9cd7a4', '#97d5a4', '#8fd2a4', '#86cfa5', '#81cda5', '#79c9a5', '#74c7a5', '#6bc4a5', '#66c2a5', '#60bba8', '#5cb7aa', '#56b0ad', '#50a9af', '#4ba4b1', '#459eb4', '#4199b6', '#3b92b9', '#378ebb', '#3387bc', '#3682ba', '#3b7cb7', '#4175b4', '#4471b2', '#496aaf', '#4d65ad', '#525fa9', '#555aa7', '#5b53a4']
color_list.extend(append_color)



## 重复

repeat_num = 30
repeat_unit = copy.deepcopy(color_list)
for i in range(repeat_num):
    color_list.extend(repeat_unit)