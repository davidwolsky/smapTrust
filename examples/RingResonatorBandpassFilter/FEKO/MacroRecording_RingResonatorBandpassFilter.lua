--CADFEKO14.0-260479 (x64)
app = cf.GetApplication()
project = app.Project

-- New project
project = app:NewProject()

-- Save project
-- app:SaveAs([[/home/rib/Documents/masters/smap/smap/examples/RingResonatorBandpassFilter/FEKO/RingResonatorBandpassFilter.cfx]])

-- Modified solution entity: Model unit
properties = project:GetProperties()
properties.ModelAttributes.Unit = cf.Enums.ModelUnitEnum.Millimetres
properties.ModelUnit = cf.Enums.ModelUnitEnum.Millimetres
project:SetProperties(properties)

-- Added variable "L1" = 24.74
L1 = project.Variables:Add("L1", "24.74")

-- Added variable "L2" = 19.51
L2 = project.Variables:Add("L2", "19.51")

-- Added variable "L3" = 24.10
L3 = project.Variables:Add("L3", "24.10")

-- Added variable "S1" = 0.293
S1 = project.Variables:Add("S1", "0.293")

-- Added variable "S2" = 0.173
S2 = project.Variables:Add("S2", "0.173")

-- Added variable "W1" = 1.232
W1 = project.Variables:Add("W1", "1.232")

-- Added variable "W2" = 0.802
W2 = project.Variables:Add("W2", "0.802")

View3D = app.Views["3D view 1"]
View3D:SetViewDirection(cf.Enums.ViewDirectionEnum.Top)

View3D:ZoomToExtents()

-- Created geometry: rectangle "topBase"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "L2"
properties.Label = "topBase"
properties.Origin.U = "-L3/2"
properties.Width = "L3"
topBase = project.Geometry:AddRectangle(properties)

-- Created geometry: rectangle "topSubtraction"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "L2-2*W2"
properties.Label = "topSubtraction"
properties.Origin.U = "-(L3/2-W1)"
properties.Origin.V = "W2"
properties.Width = "L3-2*W1"
topSubtraction = project.Geometry:AddRectangle(properties)

-- Created geometry: rectangle "sideArm"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "L1+L2"
properties.Label = "sideArm"
properties.Origin.U = "L3/2+S1"
properties.Width = "W1"
sideArm = project.Geometry:AddRectangle(properties)

-- Created geometry: rectangle "connectionArm"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "L1"
properties.Label = "connectionArm"
properties.Origin.U = "L3/2+S1+W1+S2"
properties.Origin.V = "L2"
properties.Width = "W1"
connectionArm = project.Geometry:AddRectangle(properties)

-- Added variable "A" = S2+S1
A = project.Variables:Add("A", "S2+S1")

-- Modified variable "A" from (S2+S1) to (W2+W1)
A.Expression = "W2+W1"

-- Created geometry: rectangle "connectionLine"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "A"
properties.Label = "connectionLine"
properties.Origin.U = "L3/2+S1+2*W1+S2"
properties.Origin.V = "L1+L2-A"
properties.Width = "L2/2"
connectionLine = project.Geometry:AddRectangle(properties)

-- Created geometry: subtract ""
targets = { topSubtraction }
project.Geometry:Subtract(topBase, targets)

-- Created geometry: union ""
targets = { connectionArm, connectionLine }
project.Geometry:Union(targets)

-- Created geometry: union "Union1"
-- Add a copy and mirror transform
Union1 = project.Geometry["Union1"]
Union1:Duplicate()
properties = cf.Mirror.GetDefaultProperties()
properties.Plane = cf.Enums.MirrorPlaneEnum.VN

Union1:CopyAndMirror(properties)

-- Created geometry: rectangle "sideArm"
-- Add a copy and mirror transform
sideArm:Duplicate()
properties = cf.Mirror.GetDefaultProperties()
properties.Plane = cf.Enums.MirrorPlaneEnum.VN

sideArm:CopyAndMirror(properties)

-- Save project
-- app:Save()

--CADFEKO14.0-260479 (x64)
app = cf.GetApplication()
project = app.Project

-- Added variable "eps_r" = 2.1
eps_r = project.Variables:Add("eps_r", "2.1")

-- Added medium "Dielectric1"
properties = cf.Dielectric.GetDefaultProperties()
properties.Colour = "#27DF07"
properties.DielectricModelling.RelativePermittivity = "eps_r"
properties.Label = "Dielectric1"
Dielectric1 = project.Media:AddDielectric(properties)

-- Modified medium "Dielectric1" 
properties = Dielectric1:GetProperties()
properties.Label = "Substrate"
Dielectric1:SetProperties(properties)

-- Added variable "h" = 1.5e-3
h = project.Variables:Add("h", "1.5e-3")

-- The ground plane has been set to multilayer substrate.
Environment1 = project.GroundPlane
properties = Environment1:GetProperties()
properties.DefinitionMethod = cf.Enums.GroundPlaneDefinitionMethodEnum.MultilayerSubstrate
properties.Layers[2].Medium = Dielectric1
properties.Layers[2].Thickness = "h"
Environment1:SetProperties(properties)

-- Modified variable "h" from (1.5e-3) to (1.5)
h.Expression = "1.5"

View3D = app.Views["3D view 1"]
View3D:SetViewDirection(cf.Enums.ViewDirectionEnum.Bottom)

View3D:SetViewDirection(cf.Enums.ViewDirectionEnum.Top)

-- Created port "PortInput"
properties = cf.MicrostripPort.GetDefaultProperties()
Union1_1 = project.Geometry["Union1_1"]
Edge18 = Union1_1.Edges["Edge18"]
properties.Edges = {Edge18}
properties.Label = "PortInput"
PortInput = project.Ports:AddMicrostripPort(properties)

-- Created port "PortOutput"
properties = cf.MicrostripPort.GetDefaultProperties()
Union1_2 = project.Geometry["Union1_2"]
Edge18_1 = Union1_2.Edges["Edge18"]
properties.Edges = {Edge18_1}
properties.Label = "PortOutput"
PortOutput = project.Ports:AddMicrostripPort(properties)

-- Created solution entity: Load1
properties = cf.Load.GetDefaultProperties()
properties.ImpedanceReal = "50"
properties.Label = "Load1"
properties.PortTerminal = project.Ports["PortOutput"].Terminal
Load1 = project.SolutionConfigurations["StandardConfiguration1"].Loads:AddLoad(properties)

-- Created solution entity: SParameter1
SParameterConfiguration1 = project.SolutionConfigurations:AddMultiportSParameter({project.Ports["PortInput"].Terminal})
properties = cf.SParameter.GetDefaultProperties()
properties.Label = "SParameter1"
properties.PortProperties[1].PortTerminal = project.Ports["PortInput"].Terminal
SParameter1 = SParameterConfiguration1.SParameter
SParameter1:SetProperties(properties)

-- Added variable "f_min" = 1e9
f_min = project.Variables:Add("f_min", "1e9")

-- Added variable "f_max" = 3e9
f_max = project.Variables:Add("f_max", "3e9")

-- Renamed variable "f_max" to "fmax"
f_max.Name = "fmax"

-- Renamed variable "f_min" to "fmin"
f_min.Name = "fmin"

-- Added variable "Nf" = 21
Nf = project.Variables:Add("Nf", "21")

-- Set the frequency range to discrete linearly spaced frequency samples
StandardConfiguration1 = project.SolutionConfigurations["StandardConfiguration1"]
FrequencyRange1 = StandardConfiguration1.Frequency
properties = FrequencyRange1:GetProperties()
properties.End = "fmax"
properties.NumberOfDiscreteValues = "Nf"
properties.RangeType = cf.Enums.FrequencyRangeTypeEnum.LinearSpacedDiscrete
properties.Start = "fmin"
FrequencyRange1:SetProperties(properties)

View3D:ZoomToExtents()

View3D:SetViewDirection(cf.Enums.ViewDirectionEnum.Top)

View3D:SetViewDirection(cf.Enums.ViewDirectionEnum.Bottom)

View3D:SetViewDirection(cf.Enums.ViewDirectionEnum.Top)

View3D:SetViewDirection(cf.Enums.ViewDirectionEnum.Top)

-- Save project
-- app:Save()

