-- CADFEKO v2017.2.5-317065 (x64)
app = cf.GetApplication()
project = app.Project

-- Save project
app:SaveAs([[C:\Users\ddv\Dropbox\Work\Nagraads\David Wolsky\20181010_DoubleFoldedStub\DoubleFoldedStub.cfx]])

-- Added variable "fmin" = 9.5e9
fmin = project.Variables:Add("fmin", "9.5e9")

-- Added variable "fmax" = 20e9
fmax = project.Variables:Add("fmax", "20e9")

-- Modified variable "fmin" from (9.5e9) to (5e9)
fmin.Expression = "5e9"

-- Added variable "mil" = 2.54e-5
mil = project.Variables:Add("mil", "2.54e-5")

-- Added variable "w1" = 4.8*mil
w1 = project.Variables:Add("w1", "4.8*mil")

-- Added variable "w2" = 4.8*mil
w2 = project.Variables:Add("w2", "4.8*mil")

-- Added variable "l1" = 66*mil
l1 = project.Variables:Add("l1", "66*mil")

-- Added variable "l1_mil" = 66
l1_mil = project.Variables:Add("l1_mil", "66")

-- Added variable "l2_mil" = 60
l2_mil = project.Variables:Add("l2_mil", "60")

-- Added variable "l2" = l2_mil*mil
l2 = project.Variables:Add("l2", "l2_mil*mil")

-- Added variable "s_mil" = 9.5
s_mil = project.Variables:Add("s_mil", "9.5")

-- Added variable "s" = s_mil*mil
s = project.Variables:Add("s", "s_mil*mil")

-- Modified variable "l1" from (66*mil) to (l1_mil*mil)
l1.Expression = "l1_mil*mil"

-- Added variable "t_sub" = 5*mil
t_sub = project.Variables:Add("t_sub", "5*mil")

-- Added variable "eps_r" = 9.9
eps_r = project.Variables:Add("eps_r", "9.9")

-- Added medium "Substrate"
properties = cf.Dielectric.GetDefaultProperties()
properties.Colour = "#CF474F"
properties.DielectricModelling.RelativePermittivity = "eps_r"
properties.Label = "Substrate"
Substrate = project.Media:AddDielectric(properties)

-- The ground plane has been set to multilayer substrate.
Environment1 = project.GroundPlane
properties = Environment1:GetProperties()
properties.DefinitionMethod = cf.Enums.GroundPlaneDefinitionMethodEnum.MultilayerSubstrate
FreeSpace = project.Media:GetFreeSpace()
properties.Layers[1].Medium = FreeSpace
properties.Layers[2].Medium = Substrate
properties.Layers[2].Thickness = "t_sub"
properties.Layers[3].Medium = FreeSpace
properties.ZValue = "0"
Environment1:SetProperties(properties)

-- Added variable "l3" = 30*mil
l3 = project.Variables:Add("l3", "30*mil")

-- Save project
app:Save()

-- Created geometry: rectangle "Rectangle1"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "w1"
properties.Label = "Rectangle1"
properties.Origin.N = "t_sub"
properties.Origin.U = "0"
properties.Origin.V = "-w1/2"
properties.Width = "l1/2+w2+l3"
Rectangle1 = project.Geometry:AddRectangle(properties)

-- Created geometry: rectangle "Rectangle2"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "s"
properties.Label = "Rectangle2"
properties.LocalWorkplane.Origin.X = "l1/2"
properties.LocalWorkplane.Origin.Y = "w1/2"
properties.LocalWorkplane.Origin.Z = "t_sub"
properties.Width = "w2"
Rectangle2 = project.Geometry:AddRectangle(properties)

-- Modified geometry: rectangle "Rectangle2"
properties = Rectangle2:GetProperties()
properties.Depth = "s+w2"
Rectangle2:SetProperties(properties)

-- Save project
app:Save()

-- Created geometry: rectangle "Rectangle3"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "w2"
properties.Label = "Rectangle3"
properties.LocalWorkplane.Origin.X = "l1/2"
properties.LocalWorkplane.Origin.Y = "w1/2+s"
properties.LocalWorkplane.Origin.Z = "t_sub"
properties.Origin.U = "-l2"
properties.Width = "l2"
Rectangle3 = project.Geometry:AddRectangle(properties)

-- Created geometry: union "Union1"
targets = { Rectangle1, Rectangle2, Rectangle3 }
project.Geometry:Union(targets)

properties = cf.Mirror.GetDefaultProperties()
properties.Origin.N = "t_sub"
properties.Plane = cf.Enums.MirrorPlaneEnum.VN
properties.RotationN = "0.0"

-- Add a copy and mirror transform
Union1 = project.Geometry["Union1"]
Union1:CopyAndMirror(properties)

-- Add mirror transform
properties = cf.Mirror.GetDefaultProperties()
properties.Plane = cf.Enums.MirrorPlaneEnum.UN
properties.RotationN = "0.0"
Union1_1 = project.Geometry["Union1_1"]
Mirror2 = Union1_1.Transforms:AddMirror(properties)

-- Created geometry: union "Union1"
targets = { Union1_1, Union1 }
project.Geometry:Union(targets)

-- Save project
app:Save()

-- The ground plane has been set to multilayer substrate.
properties = Environment1:GetProperties()
properties.ZValue = "t_sub"
Environment1:SetProperties(properties)

-- Save project
app:Save()

-- Created port "Port1"
properties = cf.MicrostripPort.GetDefaultProperties()
Union1_2 = project.Geometry["Union1"]
Edge2_1 = Union1_2.Edges["Edge2_1"]
properties.Edges = {Edge2_1}
properties.Label = "Port1"
Port1 = project.Ports:AddMicrostripPort(properties)

-- Created port "Port2"
properties = cf.MicrostripPort.GetDefaultProperties()
Edge2_2 = Union1_2.Edges["Edge2_2"]
properties.Edges = {Edge2_2}
properties.Label = "Port2"
Port2 = project.Ports:AddMicrostripPort(properties)

-- Save project
app:Save()

-- Set the frequency range to continuous(interpolated) frequencies.
StandardConfiguration1 = project.SolutionConfigurations["StandardConfiguration1"]
FrequencyRange1 = StandardConfiguration1.Frequency
properties = FrequencyRange1:GetProperties()
properties.End = "fmax"
properties.RangeType = cf.Enums.FrequencyRangeTypeEnum.Continuous
properties.Start = "fmin"
FrequencyRange1:SetProperties(properties)

-- Created solution entity: SParameter1
SParameterConfiguration1 = project.SolutionConfigurations:AddMultiportSParameter({
	project.Ports["Port1"].Terminal, 
	project.Ports["Port2"].Terminal})
properties = cf.SParameter.GetDefaultProperties()
properties.Label = "SParameter1"
properties.PortProperties[1].Terminal = project.Ports["Port1"].Terminal
properties.PortProperties[2] = {}
properties.PortProperties[2].Active = true
properties.PortProperties[2].Impedance = "50"
properties.PortProperties[2].Terminal = project.Ports["Port2"].Terminal
properties.TouchstoneExportEnabled = true
SParameter1 = SParameterConfiguration1.SParameter
SParameter1:SetProperties(properties)
