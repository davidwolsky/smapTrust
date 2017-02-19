-- CADFEKO v14.0.432-293039 (x64)
app = cf.GetApplication()
project = app.Project

-- Save project
-- app:SaveAs([[C:\Users\19718330\Documents\GitHub\masters\ws_smap\smapTrust\examples\BandpassFilter\bandpassFilter.cfx]])

-- Added variable "freq" = 5e9
freq = project.Variables:Add("freq", "5e9")

-- Modified solution entity: Model unit
properties = project:GetProperties()
properties.ModelAttributes.Unit = cf.Enums.ModelUnitEnum.Millimetres
project:SetProperties(properties)

-- Added variable "L1" = 6.795
L1 = project.Variables:Add("L1", "6.795")

-- Added variable "L2" = 4.848
L2 = project.Variables:Add("L2", "4.848")

-- Added variable "L3" = 6.667
L3 = project.Variables:Add("L3", "6.667")

-- Added variable "L4" = 4.956
L4 = project.Variables:Add("L4", "4.956")

-- Added variable "d" = 0.6
d = project.Variables:Add("d", "0.6")

-- Added variable "D1" = 0.6
D1 = project.Variables:Add("D1", "0.6")

-- Added variable "D2" = 5
D2 = project.Variables:Add("D2", "5")

d:Delete()

-- Added variable "h" = 0.66
h = project.Variables:Add("h", "0.66")

-- Added variable "eps_r" = 9
eps_r = project.Variables:Add("eps_r", "9")

-- Renamed variable "D1" to "d1"
D1.Name = "d1"

-- Renamed variable "D2" to "d2"
D2.Name = "d2"

-- Added variable "lambda0" = c0/freq
lambda0 = project.Variables:Add("lambda0", "c0/freq")

-- Added medium "Substrate"
properties = cf.Dielectric.GetDefaultProperties()
properties.Colour = "#07BFA7"
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
properties.Layers[2].Thickness = "h"
properties.Layers[3].Medium = FreeSpace
properties.ZValue = "0"
Environment1:SetProperties(properties)

View3D = app.Views["3D view 1"]
View3D:ZoomToExtents()

-- Added variable "g" = 0.09
g = project.Variables:Add("g", "0.09")

View3D:ZoomToExtents()

-- Created geometry: rectangle "stubLeft"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "d1"
properties.Label = "stubLeft"
properties.Origin.U = "g/2"
properties.Origin.V = "-d1/2"
properties.Width = "L1+L2"
stubLeft = project.Geometry:AddRectangle(properties)

View3D:ZoomToExtents()

-- Created geometry: rectangle "stubRight"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "d1"
properties.Label = "stubRight"
properties.Origin.U = "-g/2"
properties.Origin.V = "-d1/2"
properties.Width = "-(L1+L2)"
stubRight = project.Geometry:AddRectangle(properties)

-- Modified geometry: rectangle "stubRight"
properties = stubRight:GetProperties()
properties.Width = "-(L3+L4)"
stubRight:SetProperties(properties)

-- Created geometry: rectangle "feedRight"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "-d2"
properties.Label = "feedRight"
properties.Origin.U = "-g/2-L3"
properties.Origin.V = "-d1/2"
properties.Width = "d1"
feedRight = project.Geometry:AddRectangle(properties)

-- Created geometry: rectangle "feedLeft"
properties = cf.Rectangle.GetDefaultProperties()
properties.Depth = "d2"
properties.Label = "feedLeft"
properties.Origin.U = "g/2+L2"
properties.Origin.V = "d1/2"
properties.Width = "d1"
feedLeft = project.Geometry:AddRectangle(properties)

-- Modified geometry: rectangle "feedLeft"
properties = feedLeft:GetProperties()
properties.Origin.U = "g/2+L2-d1/2"
feedLeft:SetProperties(properties)

-- Modified geometry: rectangle "feedRight"
properties = feedRight:GetProperties()
properties.Origin.U = "-g/2-L3+d1/2"
feedRight:SetProperties(properties)

-- Modified variable "L2" from (4.848) to (4.828)
L2.Expression = "4.828"

-- Modified geometry: rectangle "feedRight"
properties = feedRight:GetProperties()
properties.Origin.U = "-g/2-L3-d1/2"
feedRight:SetProperties(properties)

-- Created geometry: union "Union1"
targets = { stubLeft, feedLeft }
project.Geometry:Union(targets)

-- Created geometry: union "Union2"
targets = { stubRight, feedRight }
project.Geometry:Union(targets)

View3D:ZoomToExtents()

View3D:SetViewDirection(cf.Enums.ViewDirectionEnum.Isometric)

View3D:SetViewDirection(cf.Enums.ViewDirectionEnum.Front)

View3D:SetViewDirection(cf.Enums.ViewDirectionEnum.Top)

View3D:ZoomToExtents()

-- Created port "PortLeft"
properties = cf.MicrostripPort.GetDefaultProperties()
Union2 = project.Geometry["Union2"]
Edge11 = Union2.Edges["Edge11"]
properties.Edges = {Edge11}
properties.Label = "PortLeft"
PortLeft = project.Ports:AddMicrostripPort(properties)

-- Created port "PortRight"
properties = cf.MicrostripPort.GetDefaultProperties()
Union1 = project.Geometry["Union1"]
Edge13 = Union1.Edges["Edge13"]
properties.Edges = {Edge13}
properties.Label = "PortRight"
PortRight = project.Ports:AddMicrostripPort(properties)

-- Created solution entity: SParameter1
SParameterConfiguration1 = project.SolutionConfigurations:AddMultiportSParameter({
	project.Ports["PortLeft"].Terminal, 
	project.Ports["PortRight"].Terminal})
properties = cf.SParameter.GetDefaultProperties()
properties.Label = "SParameter1"
properties.PortProperties[1].Terminal = project.Ports["PortLeft"].Terminal
properties.PortProperties[2] = {}
properties.PortProperties[2].Active = true
properties.PortProperties[2].Impedance = "50"
properties.PortProperties[2].Terminal = project.Ports["PortRight"].Terminal
SParameter1 = SParameterConfiguration1.SParameter
SParameter1:SetProperties(properties)

-- Added variable "freq_min" = 2e9
freq_min = project.Variables:Add("freq_min", "2e9")

-- Added variable "freq_max" = 8e9
freq_max = project.Variables:Add("freq_max", "8e9")

-- Added variable "Nf" = 100
Nf = project.Variables:Add("Nf", "100")

-- Modified variable "Nf" from (100) to (101)
Nf.Expression = "101"

-- Set the frequency range to discrete linearly spaced frequency samples
FrequencyRange1 = SParameterConfiguration1.Frequency
properties = FrequencyRange1:GetProperties()
properties.End = "freq_max"
properties.NumberOfDiscreteValues = "Nf"
properties.RangeType = cf.Enums.FrequencyRangeTypeEnum.LinearSpacedDiscrete
properties.Start = "freq_min"
FrequencyRange1:SetProperties(properties)

-- Mesh the model
project.Mesher:Mesh()

-- Updating mesh parameters
MeshSettings = project.Mesher.Settings
properties = MeshSettings:GetProperties()
properties.Advanced.MinElementSize = 37.0511713132585
properties.Advanced.RefinementFactor = 62.4196350581785
properties.MeshSizeOption = cf.Enums.MeshSizeOptionEnum.Fine
MeshSettings:SetProperties(properties)

-- Mesh the model
project.Mesher:Mesh()

-- Save project
app:Save()
