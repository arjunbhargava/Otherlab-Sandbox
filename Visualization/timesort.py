#!/usr/bin/env python

motion1 = [0, .25, .5, .75, 1]
motion2 = [0, .4, .45, .9, 1]
temp_motion1 = []
temp_motion2 = []
pos1 = 0
pos2 = 0
flag = 0

while flag == 0:
  temp_motion1.append(motion1[pos1])
  temp_motion2.append(motion2[pos2])
  if pos1 < len(motion1)-1:
    print "yeah"
    pos1 = pos1 + 1
  if pos2 < len(motion2)-1:
    print "buddy"
    pos2 = pos2 + 1
  if temp_motion1[-1]==1 and temp_motion2[-1]==1:
    flag = 1
  if pos1 <= len(motion1) - 1 and pos2 <= len(motion2) -1:
    if motion1[pos1] > motion2[pos2]:
      pos1 = pos1 - 1
    elif motion2[pos2] > motion1[pos1]: 
      pos2 = pos2 - 1