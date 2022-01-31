#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT

"""
import numpy as np
nc=64
Number = 3
Red    = np.array([ 0.00, 1.00, 0.80], np.float64)
Blue   = np.array([ 0.80, 1.00, 0.00], np.float64)
Green  = np.array([ 0.00, 1.00, 0.00], np.float64)
Length = np.array([ 0.00, 0.50, 1.00], np.float64)
FI = ROOT.TColor().CreateGradientColorTable(Number,Length,Red,Green,Blue,nc)
#MyPalette = np.fromiter([FI+i for i in range(nc)], 'int32')
#ROOT.gStyle.SetPalette(ROOT.kDarkBodyRadiator)
"""

ROOT.gStyle.SetOptStat('ne')
#ROOT.gROOT.SetStyle("Plain")
#ROOT.gROOT.ForceStyle()

ROOT.gStyle.SetPalette(56)
#ROOT.TColor.InvertPalette()


rebin_x = 50
rebin_y = 4
files = ['HorPos99.root', 'HorNeg99.root']
HDe, HDp = [], []

for i in range(2):
  XYD   = ROOT.TList()
  hfile = ROOT.TFile.Open(files[i], 'READ')
  name  = hfile.GetListOfKeys()[0].GetName()
  XYD   = hfile.Get(name)
  hfile.Close()
  HDp.append(XYD[0])
  HDe.append(XYD[1])
  HDe[i].RebinX(rebin_x)
  HDp[i].RebinX(rebin_x)
  HDe[i].RebinY(rebin_y)
  HDp[i].RebinY(rebin_y)


HDp[0].Add(HDp[1], -1)
HDe[0].Add(HDe[1], -1)
YPp = HDp[0].ProjectionX()
YPe = HDe[0].ProjectionX()

cv = ROOT.TCanvas('cv','cv',10,10,1200,1000)
cv.Divide(2,2)

HDp[0].GetXaxis().SetLabelSize(0.05); HDp[0].GetYaxis().SetLabelSize(0.05)
HDe[0].GetXaxis().SetLabelSize(0.05);
HDe[0].GetYaxis().SetLabelSize(0.05); HDe[0].GetYaxis().SetDecimals()
YPp.GetXaxis().SetLabelSize(   0.05); YPp.GetXaxis().SetDecimals()
YPp.GetYaxis().SetLabelSize(   0.05); YPp.SetLineWidth(3)
YPe.GetXaxis().SetLabelSize(   0.05); YPe.GetYaxis().SetLabelSize(0.05)
YPe.SetLineWidth(3)

cv.cd(1); cv.GetPad(1).SetGrid(); HDp[0].Draw('COLZ')
cv.cd(2); cv.GetPad(2).SetGrid(); HDe[0].Draw('COLZ')
cv.cd(3); cv.GetPad(3).SetGrid(); YPp.Draw('HIST')
cv.cd(4); cv.GetPad(4).SetGrid(); YPe.Draw('HIST')

"""
cv.cd(1); cv.GetPad(1).SetGrid(); HDp[0].Draw('COLZ')
HDp[0].GetXaxis().SetLabelSize(0.05); HDp[0].GetYaxis().SetLabelSize(0.05)
cv.cd(2); cv.GetPad(2).SetGrid(); HDe[0].Draw('COLZ')
HDe[0].GetXaxis().SetLabelSize(0.05);
HDe[0].GetYaxis().SetLabelSize(0.05); HDe[0].GetYaxis().SetDecimals()
cv.cd(3); cv.GetPad(3).SetGrid(); YPp.Draw('HIST'); YPp.SetLineWidth(3)
YPp.GetXaxis().SetLabelSize(0.05); YPp.GetXaxis().SetDecimals()
YPp.GetYaxis().SetLabelSize(0.05)
cv.cd(4); cv.GetPad(4).SetGrid(); YPe.Draw('HIST'); YPe.SetLineWidth(3)
YPe.GetXaxis().SetLabelSize(0.05); YPe.GetYaxis().SetLabelSize(0.05)
"""
cv.Modified()
cv.Update()

input()

exit()

