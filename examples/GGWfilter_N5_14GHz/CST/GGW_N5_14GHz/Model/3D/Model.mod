'# MWS Version: Version 2018.4 - Apr 06 2018 - ACIS 27.0.2 -

'# length = mm
'# frequency = GHz
'# time = ns
'# frequency range: fmin = 13.8 fmax = 14.2
'# created = '[VERSION]2017.5|26.0.1|20170804[/VERSION]


'@ define units

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Units 
     .Geometry "mm" 
     .Frequency "GHz" 
     .Time "ns" 
     .TemperatureUnit "Kelvin" 
     .Voltage "V" 
     .Current "A" 
     .Resistance "Ohm" 
     .Conductance "Siemens" 
     .Capacitance "PikoF" 
     .Inductance "NanoH" 
End With

'@ define frequency range

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solver.FrequencyRange "13.7", "14.3"

'@ activate local coordinates

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
WCS.ActivateWCS "local"

'@ move wcs

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
WCS.MoveWCS "local", "-2*W-4*a-S12-S23", "0.0", "0.0"

'@ new component: component1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Component.New "component1"

'@ define brick: component1:PinSide1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Brick
     .Reset 
     .Name "PinSide1" 
     .Component "component1" 
     .Material "PEC" 
     .Xrange "-W/2-a", "-W/2" 
     .Yrange "-a/2", "a/2" 
     .Zrange "0", "d" 
     .Create
End With

'@ transform: translate component1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1" 
     .Vector "0", "p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "2" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: mirror component1:PinSide1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinSide1_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:PinSide1_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinSide1_2" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ boolean add shapes: component1:PinSide1, component1:PinSide1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinSide1", "component1:PinSide1_1"

'@ boolean add shapes: component1:PinSide1, component1:PinSide1_1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinSide1", "component1:PinSide1_1_1"

'@ boolean add shapes: component1:PinSide1, component1:PinSide1_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinSide1", "component1:PinSide1_2"

'@ boolean add shapes: component1:PinSide1, component1:PinSide1_2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinSide1", "component1:PinSide1_2_1"

'@ transform: mirror component1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: translate component1:PinSide1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinSide1" 
     .Vector "-p-a", "0", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ define brick: component1:PinCorner1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Brick
     .Reset 
     .Name "PinCorner1" 
     .Component "component1" 
     .Material "PEC" 
     .Xrange "-W/2-a", "-W/2" 
     .Yrange "L1/2", "L1/2+a" 
     .Zrange "0", "d" 
     .Create
End With

'@ transform: translate component1:PinCorner1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinCorner1" 
     .Vector "-p-a", "0", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: translate component1:PinCorner1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinCorner1" 
     .Vector "0", "p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: translate component1:PinCorner1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinCorner1_1" 
     .Vector "0", "p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: mirror component1:PinCorner1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinCorner1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:PinCorner1_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinCorner1_2" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ define brick: component1:PinTop1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Brick
     .Reset 
     .Name "PinTop1" 
     .Component "component1" 
     .Material "PEC" 
     .Xrange "-a/2", "a/2" 
     .Yrange "L1/2", "L1/2+a" 
     .Zrange "0", "d" 
     .Create
End With

'@ transform: translate component1:PinTop1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop1" 
     .Vector "p+a", "0", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "2" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: mirror component1:PinTop1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop1_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:PinTop1_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop1_2" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: translate component1:PinTop1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop1" 
     .Vector "0", "p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: translate component1:PinTop1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop1_1" 
     .Vector "0", "p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: translate component1:PinTop1_1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop1_1_1" 
     .Vector "0", "p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: translate component1:PinTop1_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop1_2" 
     .Vector "0", "p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: translate component1:PinTop1_2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop1_2_1" 
     .Vector "0", "p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ boolean add shapes: component1:PinTop1, component1:PinTop1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop1", "component1:PinTop1_1"

'@ boolean add shapes: component1:PinTop1, component1:PinTop1_1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop1", "component1:PinTop1_1_1"

'@ boolean add shapes: component1:PinTop1, component1:PinTop1_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop1", "component1:PinTop1_2"

'@ boolean add shapes: component1:PinTop1, component1:PinTop1_2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop1", "component1:PinTop1_2_1"

'@ boolean add shapes: component1:PinTop1_1_1_1, component1:PinTop1_1_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop1_1_1_1", "component1:PinTop1_1_2"

'@ boolean add shapes: component1:PinTop1_1_1_1, component1:PinTop1_2_1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop1_1_1_1", "component1:PinTop1_2_1_1"

'@ boolean add shapes: component1:PinTop1_1_1_1, component1:PinTop1_2_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop1_1_1_1", "component1:PinTop1_2_2"

'@ boolean add shapes: component1:PinTop1_1_1_1, component1:PinTop1_3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop1_1_1_1", "component1:PinTop1_3"

'@ boolean add shapes: component1:PinCorner1, component1:PinCorner1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner1", "component1:PinCorner1_1"

'@ boolean add shapes: component1:PinCorner1, component1:PinCorner1_1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner1", "component1:PinCorner1_1_1"

'@ boolean add shapes: component1:PinCorner1, component1:PinCorner1_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner1", "component1:PinCorner1_2"

'@ boolean add shapes: component1:PinCorner1_2_1, component1:PinCorner1_3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner1_2_1", "component1:PinCorner1_3"

'@ boolean add shapes: component1:PinCorner1, component1:PinCorner1_2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner1", "component1:PinCorner1_2_1"

'@ boolean add shapes: component1:PinCorner1, component1:PinTop1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner1", "component1:PinTop1"

'@ boolean add shapes: component1:PinCorner1, component1:PinTop1_1_1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner1", "component1:PinTop1_1_1_1"

'@ transform: mirror component1:PinCorner1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinCorner1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ define material: Teflon (PTFE) (loss free)

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Material
     .Reset
     .Name "Teflon (PTFE) (loss free)"
     .Folder ""
.FrqType "all"
.Type "Normal"
.SetMaterialUnit "GHz", "mm"
.Epsilon "2.1"
.Mu "1.0"
.Kappa "0.0"
.TanD "0.0"
.TanDFreq "0.0"
.TanDGiven "False"
.TanDModel "ConstTanD"
.KappaM "0.0"
.TanDM "0.0"
.TanDMFreq "0.0"
.TanDMGiven "False"
.TanDMModel "ConstKappa"
.DispModelEps "None"
.DispModelMu "None"
.DispersiveFittingSchemeEps "General 1st"
.DispersiveFittingSchemeMu "General 1st"
.UseGeneralDispersionEps "False"
.UseGeneralDispersionMu "False"
.Rho "2200.0"
.ThermalType "Normal"
.ThermalConductivity "0.2"
.HeatCapacity "1.0"
.SetActiveMaterial "all"
.MechanicsType "Isotropic"
.YoungsModulus "0.5"
.PoissonsRatio "0.4"
.ThermalExpansionRate "140"
.Colour "0.75", "0.95", "0.85"
.Wireframe "False"
.Transparency "0"
.Create
End With

'@ move wcs

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
WCS.MoveWCS "local", "0.0", "-W/2+x", "d+h"

'@ define cylinder: component1:FeedPin

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Cylinder 
     .Reset 
     .Name "FeedPin" 
     .Component "component1" 
     .Material "PEC" 
     .OuterRadius "r" 
     .InnerRadius "0.0" 
     .Axis "z" 
     .Zrange "-pd", "coax_len" 
     .Xcenter "0" 
     .Ycenter "0" 
     .Segments "0" 
     .Create 
End With

'@ define cylinder: component1:FeedDielectric

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Cylinder 
     .Reset 
     .Name "FeedDielectric" 
     .Component "component1" 
     .Material "Teflon (PTFE) (loss free)" 
     .OuterRadius "r_out" 
     .InnerRadius "r" 
     .Axis "z" 
     .Zrange "0", "coax_len" 
     .Xcenter "0" 
     .Ycenter "0" 
     .Segments "0" 
     .Create 
End With

'@ align wcs with global coordinates

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
WCS.AlignWCSWithGlobalCoordinates

'@ move wcs

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
WCS.MoveWCS "local", "-W-2*a-S23", "0", "0.0"

'@ boolean add shapes: component1:PinCorner1, component1:PinCorner1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner1", "component1:PinCorner1_1"

'@ boolean add shapes: component1:PinCorner1, component1:PinSide1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner1", "component1:PinSide1"

'@ boolean add shapes: component1:PinCorner1, component1:PinSide1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner1", "component1:PinSide1_1"

'@ boolean add shapes: component1:PinCorner1, component1:PinSide1_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner1", "component1:PinSide1_2"

'@ rename block: component1:PinCorner1 to: component1:Pin1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Rename "component1:PinCorner1", "Pin1"

'@ define brick: component1:PinSide2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Brick
     .Reset 
     .Name "PinSide2" 
     .Component "component1" 
     .Material "PEC" 
     .Xrange "-W/2-a", "-W/2" 
     .Yrange "-a/2", "a/2" 
     .Zrange "0", "d" 
     .Create
End With

'@ clear picks

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Pick.ClearAllPicks

'@ transform: translate component1:PinSide2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinSide2" 
     .Vector "0", "p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "2" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: mirror component1:PinSide2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinSide2_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:PinSide2_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinSide2_2" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ boolean add shapes: component1:PinSide2, component1:PinSide2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinSide2", "component1:PinSide2_1"

'@ boolean add shapes: component1:PinSide2, component1:PinSide2_1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinSide2", "component1:PinSide2_1_1"

'@ boolean add shapes: component1:PinSide2, component1:PinSide2_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinSide2", "component1:PinSide2_2"

'@ boolean add shapes: component1:PinSide2, component1:PinSide2_2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinSide2", "component1:PinSide2_2_1"

'@ transform: mirror component1:PinSide2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinSide2" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ boolean add shapes: component1:PinSide2, component1:PinSide2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinSide2", "component1:PinSide2_1"

'@ define brick: component1:PinCorner2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Brick
     .Reset 
     .Name "PinCorner2" 
     .Component "component1" 
     .Material "PEC" 
     .Xrange "-W/2-a", "-W/2" 
     .Yrange "L2/2", "L2/2+A" 
     .Zrange "0", "d" 
     .Create
End With

'@ define brick: component1:PinCorner2Top

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Brick
     .Reset 
     .Name "PinCorner2Top" 
     .Component "component1" 
     .Material "PEC" 
     .Xrange "-W/2-a", "-W/2" 
     .Yrange "L1/2+p+a", "L1/2+p+2*a" 
     .Zrange "0", "d" 
     .Create
End With

'@ transform: mirror component1:PinCorner2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinCorner2" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:PinCorner2Top

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinCorner2Top" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ boolean add shapes: component1:PinCorner2, component1:PinCorner2Top

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner2", "component1:PinCorner2Top"

'@ boolean add shapes: component1:PinCorner2, component1:PinCorner2Top_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner2", "component1:PinCorner2Top_1"

'@ boolean add shapes: component1:PinCorner2, component1:PinCorner2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner2", "component1:PinCorner2_1"

'@ define brick: component1:PinTop2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Brick
     .Reset 
     .Name "PinTop2" 
     .Component "component1" 
     .Material "PEC" 
     .Xrange "-a/2", "a/2" 
     .Yrange "L2/2", "L2/2+a" 
     .Zrange "0", "d" 
     .Create
End With

'@ transform: translate component1:PinTop2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop2" 
     .Vector "p+a", "0", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "2" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: mirror component1:PinTop2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop2_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:PinTop2_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop2_2" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ boolean add shapes: component1:PinTop2, component1:PinTop2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop2", "component1:PinTop2_1"

'@ boolean add shapes: component1:PinTop2, component1:PinTop2_1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop2", "component1:PinTop2_1_1"

'@ boolean add shapes: component1:PinTop2, component1:PinTop2_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop2", "component1:PinTop2_2"

'@ boolean add shapes: component1:PinTop2, component1:PinTop2_2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop2", "component1:PinTop2_2_1"

'@ transform: translate component1:PinTop2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop2" 
     .Vector "0", "-L2/2+L1/2+p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ boolean add shapes: component1:PinTop2, component1:PinTop2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop2", "component1:PinTop2_1"

'@ transform: mirror component1:PinCorner2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinCorner2" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:PinTop2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop2" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ boolean add shapes: component1:PinCorner2, component1:PinCorner2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner2", "component1:PinCorner2_1"

'@ boolean add shapes: component1:PinCorner2, component1:PinSide2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner2", "component1:PinSide2"

'@ boolean add shapes: component1:PinCorner2, component1:PinTop2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner2", "component1:PinTop2"

'@ boolean add shapes: component1:PinCorner2, component1:PinTop2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner2", "component1:PinTop2_1"

'@ rename block: component1:PinCorner2 to: component1:Pin2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Rename "component1:PinCorner2", "Pin2"

'@ align wcs with global coordinates

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
WCS.AlignWCSWithGlobalCoordinates

'@ define brick: component1:PinSide3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Brick
     .Reset 
     .Name "PinSide3" 
     .Component "component1" 
     .Material "PEC" 
     .Xrange "-W/2-a", "-W/2" 
     .Yrange "-a/2", "a/2" 
     .Zrange "0", "d" 
     .Create
End With

'@ transform: translate component1:PinSide3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinSide3" 
     .Vector "0", "p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "2" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: mirror component1:PinSide3_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinSide3_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:PinSide3_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinSide3_2" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ boolean add shapes: component1:PinSide3, component1:PinSide3_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinSide3", "component1:PinSide3_1"

'@ boolean add shapes: component1:PinSide3, component1:PinSide3_1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinSide3", "component1:PinSide3_1_1"

'@ boolean add shapes: component1:PinSide3, component1:PinSide3_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinSide3", "component1:PinSide3_2"

'@ boolean add shapes: component1:PinSide3, component1:PinSide3_2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinSide3", "component1:PinSide3_2_1"

'@ define brick: component1:PinCorner3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Brick
     .Reset 
     .Name "PinCorner3" 
     .Component "component1" 
     .Material "PEC" 
     .Xrange "-W/2-a", "-W/2" 
     .Yrange "L3/2", "L3/2+a" 
     .Zrange "0", "d" 
     .Create
End With

'@ define brick: component1:PinTop3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Brick
     .Reset 
     .Name "PinTop3" 
     .Component "component1" 
     .Material "PEC" 
     .Xrange "-a/2", "a/2" 
     .Yrange "L3/2", "L3/2+a" 
     .Zrange "0", "d" 
     .Create
End With

'@ transform: translate component1:PinTop3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop3" 
     .Vector "-p-a", "0", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "2" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: translate component1:PinCorner3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinCorner3" 
     .Vector "0", "-L3/2+L1/2+p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: translate component1:PinTop3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop3" 
     .Vector "0", "-L3/2+L1/2+p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: translate component1:PinTop3_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop3_1" 
     .Vector "0", "-L3/2+L1/2+p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ transform: translate component1:PinTop3_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop3_2" 
     .Vector "0", "-L3/2+L1/2+p+a", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Translate" 
End With

'@ boolean add shapes: component1:PinCorner3, component1:PinCorner3_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner3", "component1:PinCorner3_1"

'@ boolean add shapes: component1:PinTop3_1, component1:PinTop3_1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop3_1", "component1:PinTop3_1_1"

'@ boolean add shapes: component1:PinTop3_1, component1:PinTop3_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop3_1", "component1:PinTop3_2"

'@ boolean add shapes: component1:PinTop3_1, component1:PinTop3_2_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop3_1", "component1:PinTop3_2_1"

'@ transform: mirror component1:PinCorner3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinCorner3" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:PinTop3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop3" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:PinTop3_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop3_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:PinTop3_3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinTop3_3" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ boolean add shapes: component1:PinCorner3, component1:PinCorner3_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner3", "component1:PinCorner3_1"

'@ boolean add shapes: component1:PinCorner3, component1:PinTop3_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner3", "component1:PinTop3_1"

'@ boolean add shapes: component1:PinCorner3, component1:PinTop3_1_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner3", "component1:PinTop3_1_1"

'@ boolean add shapes: component1:PinCorner3, component1:PinSide3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner3", "component1:PinSide3"

'@ boolean add shapes: component1:PinTop3, component1:PinTop3_2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop3", "component1:PinTop3_2"

'@ boolean add shapes: component1:PinTop3, component1:PinTop3_3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop3", "component1:PinTop3_3"

'@ boolean add shapes: component1:PinTop3, component1:PinTop3_3_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinTop3", "component1:PinTop3_3_1"

'@ transform: mirror component1:Pin1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:Pin1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:Pin2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:Pin2" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:PinCorner3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:PinCorner3" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ boolean add shapes: component1:PinCorner3, component1:PinCorner3_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner3", "component1:PinCorner3_1"

'@ boolean add shapes: component1:PinCorner3, component1:PinTop3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:PinCorner3", "component1:PinTop3"

'@ rename block: component1:PinCorner3 to: component1:Pin3

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Rename "component1:PinCorner3", "Pin3"

'@ rename block: component1:Pin1_1 to: component1:Pin5

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Rename "component1:Pin1_1", "Pin5"

'@ rename block: component1:Pin2_1 to: component1:Pin4

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Rename "component1:Pin2_1", "Pin4"

'@ transform: mirror component1:FeedDielectric

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:FeedDielectric" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:FeedPin

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:FeedPin" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:FeedDielectric_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:FeedDielectric_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "False" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Transform "Shape", "Mirror" 
End With

'@ transform: mirror component1:FeedPin_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:FeedPin_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "False" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Transform "Shape", "Mirror" 
End With

'@ define brick: component1:Cavity

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Brick
     .Reset 
     .Name "Cavity" 
     .Component "component1" 
     .Material "Vacuum" 
     .Xrange "-2.5*W-S12-S23-5*a-2*(p+a)", "0" 
     .Yrange "-L1/2-a-2*(p+a)", "0" 
     .Zrange "0", "d+h" 
     .Create
End With

'@ define material colour: Vacuum

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Material 
     .Name "Vacuum"
     .Folder ""
     .Colour "0.5", "0.8", "1" 
     .Wireframe "False" 
     .Reflection "False" 
     .Allowoutline "True" 
     .Transparentoutline "False" 
     .Transparency "80" 
     .ChangeColour 
End With

'@ transform: mirror component1:Cavity

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:Cavity" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ boolean add shapes: component1:Cavity, component1:Cavity_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:Cavity", "component1:Cavity_1"

'@ transform: mirror component1:Cavity

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Transform 
     .Reset 
     .Name "component1:Cavity" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .Material "" 
     .Transform "Shape", "Mirror" 
End With

'@ boolean add shapes: component1:Cavity, component1:Cavity_1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solid.Add "component1:Cavity", "component1:Cavity_1"

'@ activate global coordinates

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
WCS.ActivateWCS "global"

'@ pick face

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Pick.PickFaceFromId "component1:FeedDielectric", "3"

'@ define port: 1

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Port 
     .Reset 
     .PortNumber "1" 
     .Label "" 
     .NumberOfModes "1" 
     .AdjustPolarization "False" 
     .PolarizationAngle "0.0" 
     .ReferencePlaneDistance "-coax_len" 
     .TextSize "50" 
     .TextMaxLimit "0" 
     .Coordinates "Picks" 
     .Orientation "positive" 
     .PortOnBound "True" 
     .ClipPickedPortToBound "False" 
     .Xrange "-39.9", "-35.8" 
     .Yrange "-6.95", "-2.85" 
     .Zrange "11.25", "11.25" 
     .XrangeAdd "0.0", "0.0" 
     .YrangeAdd "0.0", "0.0" 
     .ZrangeAdd "0.0", "0.0" 
     .SingleEnded "False" 
     .Create 
End With

'@ pick face

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Pick.PickFaceFromId "component1:FeedDielectric_1", "3"

'@ define port: 2

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With Port 
     .Reset 
     .PortNumber "2" 
     .Label "" 
     .NumberOfModes "1" 
     .AdjustPolarization "False" 
     .PolarizationAngle "0.0" 
     .ReferencePlaneDistance "-coax_len" 
     .TextSize "50" 
     .TextMaxLimit "0" 
     .Coordinates "Picks" 
     .Orientation "positive" 
     .PortOnBound "True" 
     .ClipPickedPortToBound "False" 
     .Xrange "35.8", "39.9" 
     .Yrange "2.85", "6.95" 
     .Zrange "11.25", "11.25" 
     .XrangeAdd "0.0", "0.0" 
     .YrangeAdd "0.0", "0.0" 
     .ZrangeAdd "0.0", "0.0" 
     .SingleEnded "False" 
     .Create 
End With

'@ change solver type

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
ChangeSolverType "HF Frequency Domain"

'@ define frequency domain solver parameters

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Mesh.SetCreator "High Frequency" 
With FDSolver
     .Reset 
     .SetMethod "Tetrahedral", "Fast reduced order model" 
     .OrderTet "Second" 
     .OrderSrf "First" 
     .Stimulation "All", "All" 
     .ResetExcitationList 
     .AutoNormImpedance "False" 
     .NormingImpedance "50" 
     .ModesOnly "False" 
     .ConsiderPortLossesTet "True" 
     .SetShieldAllPorts "False" 
     .AccuracyHex "1e-6" 
     .AccuracyTet "1e-4" 
     .AccuracySrf "1e-3" 
     .LimitIterations "False" 
     .MaxIterations "0" 
     .SetCalculateExcitationsInParallel "True", "False", "" 
     .StoreAllResults "False" 
     .StoreResultsInCache "False" 
     .UseHelmholtzEquation "True" 
     .LowFrequencyStabilization "True" 
     .Type "Auto" 
     .MeshAdaptionHex "False" 
     .MeshAdaptionTet "True" 
     .AcceleratedRestart "True" 
     .FreqDistAdaptMode "Distributed" 
     .NewIterativeSolver "True" 
     .TDCompatibleMaterials "False" 
     .ExtrudeOpenBC "False" 
     .SetOpenBCTypeHex "Default" 
     .SetOpenBCTypeTet "Default" 
     .AddMonitorSamples "True" 
     .CalcStatBField "False" 
     .CalcPowerLoss "True" 
     .CalcPowerLossPerComponent "False" 
     .StoreSolutionCoefficients "True" 
     .UseDoublePrecision "False" 
     .UseDoublePrecision_ML "True" 
     .MixedOrderSrf "False" 
     .MixedOrderTet "False" 
     .PreconditionerAccuracyIntEq "0.15" 
     .MLFMMAccuracy "Default" 
     .MinMLFMMBoxSize "0.20" 
     .UseCFIEForCPECIntEq "true" 
     .UseFastRCSSweepIntEq "false" 
     .UseSensitivityAnalysis "False" 
     .SetStopSweepIfCriterionMet "True" 
     .SetSweepThreshold "S-Parameters", "0.01" 
     .UseSweepThreshold "S-Parameters", "True" 
     .SetSweepThreshold "Probes", "0.05" 
     .UseSweepThreshold "Probes", "True" 
     .SweepErrorChecks "2" 
     .SweepMinimumSamples "3" 
     .SweepConsiderAll "True" 
     .SweepConsiderReset 
     .SetNumberOfResultDataSamples "1001" 
     .SetResultDataSamplingMode "Automatic" 
     .SweepWeightEvanescent "1.0" 
     .AccuracyROM "1e-4" 
     .AddSampleInterval "", "", "1", "Automatic", "True" 
     .AddSampleInterval "", "", "", "Automatic", "False" 
     .SetUseFastResonantForSweepTet "True" 
     .MPIParallelization "False"
     .UseDistributedComputing "False"
     .NetworkComputingStrategy "RunRemote"
     .NetworkComputingJobCount "3"
     .UseParallelization "True"
     .MaxCPUs "96"
     .MaximumNumberOfCPUDevices "2"
End With
With IESolver
     .Reset 
     .UseFastFrequencySweep "True" 
     .UseIEGroundPlane "False" 
     .SetRealGroundMaterialName "" 
     .CalcFarFieldInRealGround "False" 
     .RealGroundModelType "Auto" 
     .PreconditionerType "Auto" 
     .ExtendThinWireModelByWireNubs "False" 
End With
With IESolver
     .SetFMMFFCalcStopLevel "0" 
     .SetFMMFFCalcNumInterpPoints "6" 
     .UseFMMFarfieldCalc "True" 
     .SetCFIEAlpha "0.500000" 
     .LowFrequencyStabilization "False" 
     .LowFrequencyStabilizationML "True" 
     .Multilayer "False" 
     .SetiMoMACC_I "0.0001" 
     .SetiMoMACC_M "0.0001" 
     .DeembedExternalPorts "True" 
     .SetOpenBC_XY "True" 
     .OldRCSSweepDefintion "False" 
     .SetAccuracySetting "Custom" 
     .CalculateSParaforFieldsources "True" 
     .NumberOfModesCMA "3" 
     .StartFrequencyCMA "-1.0" 
     .SetAccuracySettingCMA "Default" 
     .FrequencySamplesCMA "0" 
     .SetMemSettingCMA "Auto" 
End With

'@ define frequency range

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Solver.FrequencyRange "13.8", "14.2"

'@ define frequency domain solver parameters

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Mesh.SetCreator "High Frequency" 
With FDSolver
     .Reset 
     .SetMethod "Hexahedral", "Fast reduced order model" 
     .OrderTet "Second" 
     .OrderSrf "First" 
     .Stimulation "All", "All" 
     .ResetExcitationList 
     .AutoNormImpedance "False" 
     .NormingImpedance "50" 
     .ModesOnly "False" 
     .ConsiderPortLossesTet "True" 
     .SetShieldAllPorts "False" 
     .AccuracyHex "1e-6" 
     .AccuracyTet "1e-4" 
     .AccuracySrf "1e-3" 
     .LimitIterations "False" 
     .MaxIterations "0" 
     .SetCalculateExcitationsInParallel "True", "False", "" 
     .StoreAllResults "False" 
     .StoreResultsInCache "False" 
     .UseHelmholtzEquation "True" 
     .LowFrequencyStabilization "True" 
     .Type "Auto" 
     .MeshAdaptionHex "True" 
     .MeshAdaptionTet "True" 
     .AcceleratedRestart "True" 
     .FreqDistAdaptMode "Distributed" 
     .NewIterativeSolver "True" 
     .TDCompatibleMaterials "False" 
     .ExtrudeOpenBC "False" 
     .SetOpenBCTypeHex "Default" 
     .SetOpenBCTypeTet "Default" 
     .AddMonitorSamples "True" 
     .CalcStatBField "False" 
     .CalcPowerLoss "True" 
     .CalcPowerLossPerComponent "False" 
     .StoreSolutionCoefficients "True" 
     .UseDoublePrecision "False" 
     .UseDoublePrecision_ML "True" 
     .MixedOrderSrf "False" 
     .MixedOrderTet "False" 
     .PreconditionerAccuracyIntEq "0.15" 
     .MLFMMAccuracy "Default" 
     .MinMLFMMBoxSize "0.20" 
     .UseCFIEForCPECIntEq "true" 
     .UseFastRCSSweepIntEq "false" 
     .UseSensitivityAnalysis "False" 
     .SetStopSweepIfCriterionMet "True" 
     .SetSweepThreshold "S-Parameters", "0.01" 
     .UseSweepThreshold "S-Parameters", "True" 
     .SetSweepThreshold "Probes", "0.05" 
     .UseSweepThreshold "Probes", "True" 
     .SweepErrorChecks "2" 
     .SweepMinimumSamples "3" 
     .SweepConsiderAll "True" 
     .SweepConsiderReset 
     .SetNumberOfResultDataSamples "1001" 
     .SetResultDataSamplingMode "Automatic" 
     .SweepWeightEvanescent "1.0" 
     .AccuracyROM "1e-4" 
     .AddSampleInterval "", "", "1", "Automatic", "True" 
     .AddSampleInterval "", "", "", "Automatic", "False" 
     .SetUseFastResonantForSweepTet "True" 
     .MPIParallelization "False"
     .UseDistributedComputing "False"
     .NetworkComputingStrategy "RunRemote"
     .NetworkComputingJobCount "3"
     .UseParallelization "True"
     .MaxCPUs "96"
     .MaximumNumberOfCPUDevices "2"
End With
With IESolver
     .Reset 
     .UseFastFrequencySweep "True" 
     .UseIEGroundPlane "False" 
     .SetRealGroundMaterialName "" 
     .CalcFarFieldInRealGround "False" 
     .RealGroundModelType "Auto" 
     .PreconditionerType "Auto" 
     .ExtendThinWireModelByWireNubs "False" 
End With
With IESolver
     .SetFMMFFCalcStopLevel "0" 
     .SetFMMFFCalcNumInterpPoints "6" 
     .UseFMMFarfieldCalc "True" 
     .SetCFIEAlpha "0.500000" 
     .LowFrequencyStabilization "False" 
     .LowFrequencyStabilizationML "True" 
     .Multilayer "False" 
     .SetiMoMACC_I "0.0001" 
     .SetiMoMACC_M "0.0001" 
     .DeembedExternalPorts "True" 
     .SetOpenBC_XY "True" 
     .OldRCSSweepDefintion "False" 
     .SetAccuracySetting "Custom" 
     .CalculateSParaforFieldsources "True" 
     .NumberOfModesCMA "3" 
     .StartFrequencyCMA "-1.0" 
     .SetAccuracySettingCMA "Default" 
     .FrequencySamplesCMA "0" 
     .SetMemSettingCMA "Auto" 
End With

'@ set pba mesh type

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
Mesh.MeshType "PBA"

'@ set 3d mesh adaptation results

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
With MeshSettings
   .SetMeshType "Hex"
   .Set "StepsPerWaveNear", "30"
   .Set "StepsPerWaveFar", "30"
   .Set "StepsPerBoxNear", "35"
   .Set "StepsPerBoxFar", "16"
End With

'@ deactivate high frequency solver mesh adaptation

'[VERSION]2017.5|26.0.1|20170804[/VERSION]
FDSolver.MeshAdaptionHex "False" 

'@ switch working plane

'[VERSION]2018.4|27.0.2|20180406[/VERSION]
Plot.DrawWorkplane "false" 


