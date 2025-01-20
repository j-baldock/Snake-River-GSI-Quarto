Comprehensive Aquatic Barrier Inventory - dams and assessed road-related barriers
--------------------------------------------------------

Data version: 3.14.0
Data published: 12/04/2024
Downloaded from: https://aquaticbarriers.org
Filename: aquatic_barrier_ranks.csv

Selected Area:
--------------
HUC6: 170401


Description:
------------
This inventory is a growing and living database of dams and surveyed road/stream
crossings (potential barriers) compiled by the Southeast Aquatic Resources
Partnership with the generous support from many partners and funders.
Information about network connectivity, landscape condition, and presence of
threatened and endangered aquatic organisms are added to this inventory to help
you investigate barriers at any scale for your desired purposes.

This inventory consists of datasets from local, state, and federal partners. It
is supplemented with input from partners with on the ground knowledge of
specific structures. The information on barriers is not complete or
comprehensive across the region, and depends on the availability and
completeness of existing data and level of partner feedback. Some areas of the
region are more complete than others but none should be considered 100%
complete.

All network analyses were conducted using the NHD High Resolution Plus dataset
(https://www.usgs.gov/core-science-systems/ngp/national-hydrography/nhdplus-high-resolution).

If you are able to help improve the inventory by sharing data or assisting with field
reconnaissance, please contact us (https://southeastaquatics.net/about/contact-us).

Note: data come from a variety of sources and available descriptive information
is not comprehensive.

Also note: information on rare species is highly limited.


File Contents:
--------------
BarrierType: Type of barrier
lat: Latitude in WGS84 geographic coordinates.
lon: Longitude in WGS84 geographic coordinates.
Name: dam or assessed road/stream crossing name, if available.
SARPID: SARP Identifier.
Source: Source of this record in the inventory.
SourceID: Identifier of this dam or assessed road/stream crossing in the source database
Snapped: Indicates if the dam or assessed road/stream crossing was snapped to a flowline.  Note: not all barriers snapped to flowlines are used in the network connectivity analysis.
NHDPlusID: Unique NHD Plus High Resolution flowline identifier to which this dam or assessed road/stream crossing is snapped.  -1 = not snapped to a flowline.  Note: not all barriers snapped to flowlines are used in the network connectivity analysis.
River: River or stream name where dam or assessed road/stream crossing occurs, if available.
StreamSizeClass: Stream size class based on total catchment drainage area in acres.  1a: <2,471 acres, 1b: 2,471-24,709 acres, 2: 24,710-127,999 acres, 3a: 128,000-640,001 acres, 3b: 640,002-2,471,049 acres, 4: 2,471,050-6,177,624 acres, 5: >= 6,177,625 acres.
FlowsToOcean: indicates if this dam or assessed road/stream crossing was snapped to a stream or river that is known to flow into the ocean.  Note: this underrepresents any networks that traverse regions outside the analysis region that would ultimately connect the networks to the ocean.
FlowsToGreatLakes: indicates if this dam or assessed road/stream crossing was snapped to a stream or river that is known to flow into the Great Lakes.  Note: this underrepresents any networks that traverse regions outside the analysis region that would ultimately connect the networks to the Great Lakes.
NIDID: National Inventory of Dams Identifier (legacy ID); this value was provided in earlier versions of NID or partner databases and may no longer match the latest ID, and may also have duplicate IDs or incorrectly-associated IDs.
NIDFederalID: National Inventory of Dams Federal Identifier (new ID) that can be used to join to the latest version of NID.
PartnerID: Identifier used by local partners for this dam or assessed road/stream crossing
Estimated: Dam represents an estimated dam location based on NHD high resolution waterbodies or other information.
AnnualVelocity: Annual velocity at the downstream end of the NHD Plus High Resolution flowline to which this dam or assessed road/stream crossing has been snapped, in square feet per second.  -1 if not snapped to flowline or otherwise not available
AnnualFlow: Annual flow at the downstream end of the NHD Plus High Resolution flowline to which this dam or assessed road/stream crossing has been snapped, in square cubic feet per second.  -1 if not snapped to flowline or otherwise not available
TotDASqKm: Total drainage area in square kilometers at the downstream end of the NHD Plus High Resolution flowline to which this dam or assessed road/stream crossing has been snapped, as reported by NHD.  -1 if not snapped to flowline or otherwise not available
YearCompleted: year that construction was completed, if available.  0 = data not available.
FERCRegulated: Identifies if the dam or assessed road/stream crossing is regulated by the Federal Energy Regulatory Commission, if known.
StateRegulated: Identifies if the dam or assessed road/stream crossing is regulated at the state level, if known.
NRCSDam: Identifies if the dam or assessed road/stream crossing is a flood control dam constructed in partnership with the USDA Natural Resources Conservation Service.
FedRegulatoryAgency: Identifies the federal regulatory agency for this dam or assessed road/stream crossing, if known
WaterRight: Identifies if the dam or assessed road/stream crossing has an associated water right, if known.
IsPriority: Indicates if the dam or assessed road/stream crossing has been identified as a priority by resource managers, if known.
Height: dam or assessed road/stream crossing height in feet, if available.  0 = data not available.
Length: dam or assessed road/stream crossing length in feet, if available.  0 = data not available.
Width: dam or assessed road/stream crossing width in feet, if available.  0 = data not available.
Hazard: Hazard rating of this dam or assessed road/stream crossing, if known.
Construction: material used in dam or assessed road/stream crossing construction, if known.
Purpose: primary purpose of dam or assessed road/stream crossing, if known.
Passability: passability of the dam or assessed road/stream crossing, if known.   Note: assessment dates are not known.
FishScreen: Whether or not a fish screen is present, if known.
ScreenType: Type of fish screen, if known.
Feasibility: feasibility of dam or assessed road/stream crossing removal, based on reconnaissance.  Note: reconnaissance information is available only for a small number of dam or assessed road/stream crossings.
Recon: Field reconnaissance notes, if available.
Diversion: Identifies if dam is known to be a diversion.  Note: diversion information is available only for a small number of dams.
LowheadDam: Identifies if dam is known or estimated to be a lowhead dam.  Note: lowhead dam information is available only for a small number of dams.
StorageVolume: Identifies the reported normal storage volume (acre/feet) of the impounded waterbody, from NID NormStor attribute, if known.  
WaterbodyAcres: area in acres of waterbody associated with dam.  -1 = no associated waterbody
Fatality: number of fatalities recorded at this location from Locations of Fatalities at Submerged Hydraulic Jumps (https://krcproject.groups.et.byu.net/browse.php)
Removed: Identifies if the dam or assessed road/stream crossing has been removed for conservation, if known.  Removed barriers will not have values present for all fields.
YearRemoved: year that barrier was removed or mitigated, if available.  All barriers removed prior to 2000 or where YearRemoved is unknown were lumped together for the network analysis.  0 = data not available or not removed / mitigated.
Condition: Condition of the dam or assessed road/stream crossing as of last assessment, if known. Note: assessment dates are not known.
PassageFacility: type of fish passage facility, if known.
TESpp: Number of federally-listed threatened or endangered aquatic species, compiled from element occurrence data within the same subwatershed (HUC12) as the dam or assessed road/stream crossing. Note: rare species information is based on occurrences within the same subwatershed as the barrier.  These species may or may not be impacted by this dam or assessed road/stream crossing.  Information on rare species is very limited and comprehensive information has not been provided for all states at this time.
StateSGCNSpp: Number of state-listed Species of Greatest Conservation Need (SGCN), compiled from element occurrence data within the same subwatershed (HUC12) as the dam or assessed road/stream crossing.  Note: rare species information is based on occurrences within the same subwatershed as the dam or assessed road/stream crossing.  These species may or may not be impacted by this dam or assessed road/stream crossing.  Information on rare species is very limited and comprehensive information has not been provided for all states at this time.
RegionalSGCNSpp: Number of regionally-listed Species of Greatest Conservation Need (SGCN), compiled from element occurrence data within the same subwatershed (HUC12) as the dam or assessed road/stream crossing.  Note: rare species information is based on occurrences within the same subwatershed as the dam or assessed road/stream crossing.  These species may or may not be impacted by this dam or assessed road/stream crossing.  Information on rare species is very limited and comprehensive information has not been provided for all states at this time.
Trout: Identifies one or more interior or eastern native trout species (Apache, brook, bull, cutthroat, Gila, lake, and redband) that are present within the same subwatershed (HUC12) as the dam or assessed road/stream crossing based on in available natural heritage data and other data sources.  Note: absence means that occurrences were not present in the available natural heritage data and should not be interpreted as true absences.
OwnerType: Land ownership type. This information is derived from the BLM Surface Management Agency dataset for federal lands and CBI Protected Areas Database and TNC Secured Lands Database for non-federal lands, to highlight ownership types of particular importance to partners.  NOTE: does not include most private land.
BarrierOwnerType: Barrier ownership type, if available.  For unsurveyed road / stream crossings, this information is derived from the National Bridge Inventory, US Census TIGER Roads route type, and USFS National Forest road / stream crossings database ownership information, and may not be fully accurate.
ProtectedLand: Indicates if the dam or assessed road/stream crossing occurs on public land as represented within the BLM Surface Management Agency dataset, CBI Protected Areas Database of the U.S., or TNC Secured Lands Database.
SourceLink: Link to additional information about this dam or assessed road/stream crossing provided by the data source.
EJTract: Within an overburdened and underserved Census tracts a defined by the Climate and Environmental Justice Screening tool.
EJTribal: Within a disadvantaged tribal community as defined by the Climate and Environmental Justice Screening tool based on American Indian and Alaska Native areas as defined by the US Census Bureau.  Note: all tribal communities considered disadvantaged by the Climate and Environmental Justice Screening tool.
NativeTerritories: Native / indigenous people's territories as mapped by Native Land Digital (https://native-land.ca/
WildScenicRiver: Name of the designated Wild & Scenic River(s) if within 250 meters of this dam or assessed road/stream crossing.
FishHabitatPartnership: Fish Habitat Partnerships working in the area where the dam or assessed road/stream crossing occurs.  See https://www.fishhabitat.org/the-partnerships for more information.
Basin: Name of the hydrologic basin (HUC6) where the dam or assessed road/stream crossing occurs.
Subbasin: Name of the hydrologic subbasin (HUC8) where the dam or assessed road/stream crossing occurs.
Subwatershed: Name of the hydrologic subwatershed (HUC12) where the dam or assessed road/stream crossing occurs.
HUC2: Hydrologic region identifier where the dam or assessed road/stream crossing occurs
HUC6: Hydrologic basin identifier where the dam or assessed road/stream crossing occurs.
HUC8: Hydrologic subbasin identifier where the dam or assessed road/stream crossing occurs.
HUC10: Hydrologic watershed identifier where the dam or assessed road/stream crossing occurs.
HUC12: Hydrologic subwatershed identifier where the dam or assessed road/stream crossing occurs.
State: State where dam or assessed road/stream crossing occurs.
County: County where dam or assessed road/stream crossing occurs.
CongressionalDistrict: Congressional District where dam or assessed road/stream crossing occurs (118th Congressional Districts).
Excluded: this dam or assessed road/stream crossing was excluded from the connectivity analysis based on field reconnaissance or manual review of aerial imagery.
Invasive: this dam or assessed road/stream crossing is identified as a beneficial to restricting the movement of invasive species and is not ranked
OnLoop: this dam or assessed road/stream crossing occurs on a loop within the NHD High Resolution aquatic network and is considered off-network for purposes of network analysis and ranking
HasNetwork: indicates if this dam or assessed road/stream crossing was snapped to the aquatic network for analysis.  1 = on network, 0 = off network.  Note: network metrics and scores are not available for dam or assessed road/stream crossings that are off network.
Ranked: this dam or assessed road/stream crossing was included for prioritization.  Some barriers that are beneficial to restricting the movement of invasive species or that are water diversions without associated barriers are excluded from ranking.
Intermittent: indicates if this dam or assessed road/stream crossing was snapped to a a stream or river reach coded by NHDPlusHR as an intermittent or ephemeral. -1 = not available.
Canal: indicates if this dam or assessed road/stream crossing was snapped to a a stream or river reach coded by NHDPlusHR as a canal or ditch. -1 = not available.
StreamOrder: NHDPlus Modified Strahler stream order. -1 = not available.
Landcover: average amount of the river floodplain in the upstream network that is in natural landcover types.  -1 = not available.
SizeClasses: number of unique upstream size classes that could be gained by removal of this dam or assessed road/stream crossing. -1 = not available.
PerennialSizeClasses: number of unique upstream size classes of perennial stream reaches that could be gained by removal of this dam or assessed road/stream crossing. -1 = not available.
MainstemSizeClasses: number of unique upstream size classes within mainstem networks that could be gained by removal of this dam or assessed road/stream crossing. -1 = not available.
TotalUpstreamMiles: number of miles in the upstream functional river network from this dam or assessed road/stream crossing, including miles in waterbodies. -1 = not available.
PerennialUpstreamMiles: number of perennial miles in the upstream functional river network from this dam or assessed road/stream crossing, including miles in waterbodies.  Perennial reaches are all those not specifically coded by NHD as ephemeral or intermittent, and include other types, such as canals and ditches that may not actually be perennial.  Networks are constructed using all flowlines, not just perennial reaches. -1 = not available.
IntermittentUpstreamMiles: number of ephemeral and intermittent miles in the upstream functional river network from this dam or assessed road/stream crossing, including miles in waterbodies.  Ephemeral and intermittent reaches are all those that are specifically coded by NHD as ephemeral or intermittent, and specifically excludes other types, such as canals and ditches that may actually be ephemeral or intermittent in their flow frequency.  -1 = not available.
AlteredUpstreamMiles: number of altered miles in the upstream functional river network from this dam or assessed road/stream crossing, including miles in waterbodies.  Altered reaches are those specifically identified in NHD or the National Wetlands Inventory as altered (canal / ditch, within a reservoir, or other channel alteration). -1 = not available.
UnalteredUpstreamMiles: number of unaltered miles in the upstream functional river network from this dam or assessed road/stream crossing, including miles in waterbodies.  Unaltered miles exclude reaches specifically identified in NHD or the National Wetlands Inventory as altered (canal / ditch, within a reservoir, or other channel alteration). -1 = not available.
PercentUnaltered: percent of the total upstream functional river network length from this dam or assessed road/stream crossing that is not specifically identified in NHD or the National Wetlands Inventory as altered (canal / ditch, within a reservoir, or other channel alteration).  -1 = not available.
PerennialUnalteredUpstreamMiles: number of unaltered perennial miles in the upstream functional river network from this dam or assessed road/stream crossing, including miles in waterbodies.  Unaltered miles exclude reaches specifically identified in NHD or the National Wetlands Inventory as altered (canal / ditch, within a reservoir, or other channel alteration). -1 = not available.
PercentPerennialUnaltered: percent of the perennial upstream functional river network length from this dam or assessed road/stream crossing that is not specifically identified in NHD or the National Wetlands Inventory as altered (canal / ditch, within a reservoir, or other channel alteration).  See PerennialUpstreamMiles.  -1 = not available.
PercentResilient: percent of the the upstream functional river network length from this dam or assessed road/stream crossing that is within watersheds identified by The Nature Conservancy with above average or greater freshwater resilience (v0.44), including miles in waterbodies.  -1 = not available.  See https://www.maps.tnc.org/resilientrivers/#/explore for more information.
UpstreamDrainageAcres: approximate drainage area in acres of all NHD High Resolution catchments within upstream functional network of dam or assessed road/stream crossing.  Includes the total catchment area of any NHD High Resolution flowlines that are cut by barriers in the analysis, which may overrepresent total drainage area of the network. -1 = not available.
UpstreamUnalteredWaterbodyAcres: area in acres of all unaltered lakes and ponds that intersect any reach in upstream functional network; these waterbodies are not specifically marked by data sources as altered and are not associated with dams in this inventory.  Use with caution because waterbodies may not be correctly subdivided by dams in the data sources, and altered waterbodies may not be marked as such. -1 = not available.
UpstreamUnalteredWetlandAcres: area in acres of all unaltered freshwater wetlands that intersect any reach in the upstream functional network. Wetlands are derived from specific wetland types in the National Wetlands Inventory (freshwater scrub-shrub, freshwater forested, freshwater emergent) and NHD (swamp/marsh) and exclude any specifically marked by the data provider as altered.  Use with caution because wetlands may not be correctly subdivided by dams or dikes in the data sources, and altered wetlands may not be marked as such. -1 = not available.
TotalDownstreamMiles: number of miles in the complete downstream functional river network from this dam or assessed road/stream crossing, including miles in waterbodies.  Note: this measures the length of the complete downstream network including all tributaries, and is not limited to the shortest downstream path.  -1 = not available.
FreeDownstreamMiles: number of free-flowing miles in the downstream functional river network.  Excludes miles in altered reaches in waterbodies.  -1 = not available.
FreePerennialDownstreamMiles: number of free-flowing perennial miles in the downstream functional river network.  Excludes miles in altered reaches in waterbodies.  See PerennialUpstreamMiles for definition of perennial reaches. -1 = not available.
FreeIntermittentDownstreamMiles: number of free-flowing ephemeral and intermittent miles in the downstream functional river network.  Excludes miles altered reaches in waterbodies.  See IntermittentUpstreamMiles for definition of intermittent reaches. -1 = not available.
FreeAlteredDownstreamMiles: number of free-flowing altered miles in the downstream functional river network from this dam or assessed road/stream crossing.  Excludes miles in altered reaches in waterbodies.  See AlteredUpstreaMiles for definition of altered reaches.  -1 = not available.
FreeUnalteredDownstreamMiles: number of free-flowing altered miles in the downstream functional river network from this dam or assessed road/stream crossing.  Excludes miles in altered reaches in waterbodies.  See UnalteredUpstreamMiles for definition of unaltered reaches.  -1 = not available.
GainMiles: absolute number of miles that could be gained by removal of this dam or assessed road/stream crossing.  Calculated as the minimum of the TotalUpstreamMiles and FreeDownstreamMiles unless the downstream network flows into the Great Lakes or the Ocean and has no downstream barriers, in which case this is based only on TotalUpstreamMiles. For removed barriers, this is based on the barriers present at the time this barrier was removed, with the exception of those that are immediately upstream and removed in the same year.  -1 = not available.
PerennialGainMiles: absolute number of perennial miles that could be gained by removal of this dam or assessed road/stream crossing.  Calculated as the minimum of the PerennialUpstreamMiles and FreePerennialDownstreamMiles unless the downstream network flows into the Great Lakes or the Ocean and has no downstream barriers, in which case this is based only on PerennialUpstreamMiles.  For removed barriers, this is based on the barriers present at the time this barrier was removed, with the exception of those that are immediately upstream and removed in the same year.  -1 = not available.
TotalNetworkMiles: sum of TotalUpstreamMiles and FreeDownstreamMiles. -1 = not available.
TotalPerennialNetworkMiles: sum of PerennialUpstreamMiles and FreePerennialDownstreamMiles. -1 = not available.
TotalUpstreamMainstemMiles: number of miles in the upstream mainstem river network from this dam or assessed road/stream crossing, including miles in waterbodies.  Mainstems are defined as having >= 1 square mile drainage area and on the same stream order as the reach where the dam is located. -1 = not available.
PerennialUpstreamMainstemMiles: number of perennial miles in the upstream mainstem river network from this dam or assessed road/stream crossing, including miles in waterbodies.  See TotalUpstreamMainstemMiles for definition of mainstem and PerennialUpstreamMiles for definition of perennial reaches.  Networks are constructed using all flowlines, not just perennial reaches. -1 = not available.
IntermittentUpstreamMainstemMiles: number of ephemeral and intermittent miles in the upstream mainstem river network from this dam or assessed road/stream crossing, including miles in waterbodies.  See TotalUpstreamMainstemMiles for definition of mainstem and IntermittentUpstreamMiles for definition of intermittent reaches.  -1 = not available.
AlteredUpstreamMainstemMiles: number of altered miles in the upstream mainstem river network from this dam or assessed road/stream crossing, including miles in waterbodies.  See TotalUpstreamMainstemMiles for definition of mainstem and AlteredUpstreamMiles for definition of altered reaches.  -1 = not available.
UnalteredUpstreamMainstemMiles: number of unaltered miles in the upstream mainstem river network from this dam or assessed road/stream crossing, including miles in waterbodies.  See TotalUpstreamMainstemMiles for definition of mainstem and UnalteredUpstreamMiles for definition of unaltered reaches.  -1 = not available.
PerennialUnalteredUpstreamMainstemMiles: number of unaltered perennial miles in the upstream mainstem river network from this dam or assessed road/stream crossing, including miles in waterbodies.  See TotalUpstreamMainstemMiles for definition of mainstem and PerennialUnalteredUpstreamMiles for definition of perennial unaltered reaches.   -1 = not available.
PercentMainstemUnaltered: percent of the upstream mainstem river network length from this dam or assessed road/stream crossing that is not specifically identified in NHD or the National Wetlands Inventory as altered (canal / ditch, within a reservoir, or other channel alteration).  -1 = not available.
FreeLinearDownstreamMiles: number of miles in the linear downstream flow direction between this dam or assessed road/stream crossing and the next barrier downstream (if any) or downstream-most point (e.g., ocean, river outlet, interior basin, etc) on the full aquatic network on which it occurs. Excludes miles in altered reaches in waterbodies.  -1 = not available.
FreePerennialLinearDownstreamMiles: number of perennial miles in the linear downstream flow direction between this dam or assessed road/stream crossing and the next barrier downstream or downstream-most point on the full aquatic network on which it occurs.  Excludes miles in altered reaches in waterbodies.  See PerennialUpstreamMiles for definition of perennial reaches. -1 = not available.
FreeIntermittentLinearDownstreamMiles: number of ephemeral and intermittent miles in the linear downstream flow direction between this dam or assessed road/stream crossing and the next barrier downstream or downstream-most point on the full aquatic network on which it occurs.  Excludes miles in altered reaches in waterbodies.  See IntermittentUpstreamMiles for definition of intermittent reaches. -1 = not available.
FreeAlteredLinearDownstreamMiles: number of altered miles in the linear downstream flow direction between this dam or assessed road/stream crossing and the next barrier downstream or downstream-most point on the full aquatic network on which it occurs.  Excludes miles in altered reaches in waterbodies.  See AlteredUpstreamMiles for definition of altered reaches. -1 = not available.
FreeUnalteredLinearDownstreamMiles: number of unaltered miles in the linear downstream flow direction between this dam or assessed road/stream crossing and the next barrier downstream or downstream-most point on the full aquatic network on which it occurs.  Excludes miles in unaltered reaches in waterbodies.  See UnalteredUpstreamMiles for definition of unaltered reaches. -1 = not available.
MainstemGainMiles: absolute number of mainstem miles that could be gained by removal of this dam or assessed road/stream crossing.  Calculated as the minimum of the TotalUpstreamMainstemMiles and FreeLinearDownstreamMiles unless the downstream network flows into the Great Lakes or the Ocean and has no downstream barriers, in which case this is based only on TotalUpstreamMainstemMiles.  -1 = not available.
PerennialMainstemGainMiles: absolute number of perennial mainstem miles that could be gained by removal of this dam or assessed road/stream crossing.  Calculated as the minimum of the PerennialUpstreamMainstemMiles and FreePerennialLinearDownstreamMiles unless the downstream network flows into the Great Lakes or the Ocean and has no downstream barriers, in which case this is based only on TotalUpstreamMainstemMiles.  -1 = not available.
TotalMainstemNetworkMiles: sum of TotalUpstreamMainstemMiles and FreeLinearDownstreamMiles. -1 = not available.
TotalMainstemPerennialNetworkMiles: sum of PerennialUpstreamMainstemMiles and FreePerennialLinearDownstreamMiles. -1 = not available.
UpstreamWaterfalls: number of waterfalls at the upstream ends of the functional network for this dam or assessed road/stream crossing. -1 = not available.
UpstreamDams: number of dams at the upstream ends of the functional network for this dam or assessed road/stream crossing.
UpstreamSmallBarriers: number of assessed road/stream crossings within the functional network if this barrier is a dam or at the upstream ends of the functional network if this barrier is a road/stream crossing. -1 = not available.
UpstreamRoadCrossings: number of uninventoried estimated road crossings within the functional network for this dam or assessed road/stream crossing. -1 = not available.
UpstreamHeadwaters: number of headwaters within the functional network for this dam or assessed road/stream crossing. -1 = not available.
TotalUpstreamRoadCrossings: total number of uninventroeid estimated road crossings upstream of this dam or assessed road/stream crossing; includes in all functional networks above this dam or assessed road/stream crossing. -1 = not available.
TotalDownstreamWaterfalls: total number of waterfalls between this dam or assessed road/stream crossing and the downstream-most point the full aquatic network on which it occurs. -1 = not available.
TotalDownstreamDams: total number of dams between this dam or assessed road/stream crossing and the downstream-most point the full aquatic network on which it occurs (e.g., river mouth). -1 = not available.
TotalDownstreamSmallBarriers: total number of assessed road/stream crossings between this dam or assessed road/stream crossing and the downstream-most point the full aquatic network on which it occurs. -1 = not available.
MilesToOutlet: miles between this dam or assessed road/stream crossing and the downstream-most point (e.g., ocean, river outlet, interior basin, etc) on the full aquatic network on which it occurs. -1 = not available.
InvasiveNetwork: indicates if there is an invasive species barrier at or downstream of this dam or assessed road/stream crossing.
YellowstoneCutthroatTroutHabitatUpstreamMiles: number of miles in the upstream river network from this dam or assessed road/stream crossing that are attributed as habitat for Yellowstone cutthroat trout.  Habitat reaches are not necessarily contiguous.  Habitat is estimated at the NHDPlusHR flowline level based on best available habitat data provided by StreamNet; please see https://aquaticbarriers.org/habitat_methods for more information. -1 = not available.
FreeYellowstoneCutthroatTroutHabitatDownstreamMiles: number of free-flowing miles in the downstream river network from this dam or assessed road/stream crossing that are attributed as habitat for Yellowstone cutthroat trout.  Habitat reaches are not necessarily contiguous.  Habitat is estimated at the NHDPlusHR flowline level based on best available habitat data provided by StreamNet; please see https://aquaticbarriers.org/habitat_methods for more information. -1 = not available.
CrossingCode: crossing identifier.
NearestCrossingID: The SARPID of the nearest road/stream crossing point, if any are found within 50 meters
NearestUSGSCrossingID: The ID of the nearest road/stream crossing point in the USGS Database of Stream Crossings in the United States (2022), if any are found within 50 meters.  Will be blank for crossings from other sources
Road: road name, if available.
RoadType: type of road, if available.
CrossingType: type of road / stream crossing, if known.
Constriction: type of constriction at road / stream crossing, if known.
PotentialProject: reconnaissance information about the crossing, including severity of the barrier and / or potential for removal project.
BarrierSeverity: barrier severity of the dam or assessed road/stream crossing, if known.   Note: assessment dates are not known.
SARP_Score: The best way to consider the aquatic passability scores is that they represent the degree to which crossings deviate from an ideal crossing. We assume that those crossings that are very close to the ideal (scores > 0.6) will present only a minor or insignificant barrier to aquatic organisms. Those structures that are farthest from the ideal (scores < 0.4) are likely to be either significant or severe barriers. These are, however, arbitrary distinctions imposed on a continuous scoring system and should be used with that in mind. -1 = not available.
ProtocolUsed: Name of survey protocol used