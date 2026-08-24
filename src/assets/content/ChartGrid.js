export default {
    chartGridItems: [
        //we'll be replacing these img_src with paths to location on s3
        //vizRoutes will direct to appropriate subpage
        {
            title: 'Threats to inland fisheries',
            project: 'Findex',
            vizKey: 'FindexThreatDefinitions',
            vizRoute: 'inland-fish-threats',
            img_src: 'inland_fisheries_thumbnail.webp',
            alt: '',
            chartOrder: 1,
            description: 'Inland fisheries face many threats.'
        },{
            title: 'Inland fisheries are threatened',
            project: 'Findex',
            vizKey: 'FindexThreatSankey',
            vizRoute: 'inland-fish-total-threats',
            img_src: 'findex_sankey_thumbnail.webp',
            alt: '',
            chartOrder: 2,
            description: 'Land use change is threatening inland fisheries.'
        },
        {
            title: 'Climate vulnerability of recreational fish consumption by country',
            project: 'Rec Fish as Food',
            vizKey: 'FishAsFoodClimate',
            vizRoute: 'rec-fish-vulnerability',
            img_src: 'FaF_climate_vulnerability_thumbnail.webp',
            chartOrder: 3,
            alt: 'Bubble chart showing climate vulnerability of recreational fish consumption by country, with circles colored by continent and sized by number of fishers.',
            description: 'Climate vulnerability of recreational fish consumption varies by country.'
        },
        {
            title: 'A global view of threats',
            project: 'Findex',
            vizKey: 'FindexGlobalThreats',
            vizRoute: 'inland-fish-global-threats',
            img_src: 'threat_by_basin_thumbnail.webp',
            alt: '',
            chartOrder: 3,
            description: 'Inland fisheries are threatened globally.'
        },
        {
            title: 'Explore the Juneau Icefield',
            project: 'Fire in Ice',
            vizKey: 'GlacierScan',
            vizRoute: 'glacier-scan',
            img_src: 'glacial_mri_thumbnail.webp',
            alt: 'Interactive topographic map of the Juneau Icefield showing snow core collection sites.',
            chartOrder: 1,
            description: 'Glacier ice records change.'
        },
        {
            title: 'Wildfire aerosols',
            project: 'Fire in Ice',
            vizKey: 'WildfireAerosols',
            vizRoute: 'wildfire-aerosols',
            img_src: 'wildfire_aerosols_thumbnail.webp',
            alt: 'Visual representation of a snow core showing layers of wildfire aerosol deposition preserved in glacier ice.',
            chartOrder: 2,
            description: 'Wildfire aerosols are preserved in glaciers.'
        },          
        {
            title: 'Regional wildfire deposition',
            project: 'Fire in Ice',
            vizKey: 'RegionalFires',
            vizRoute: 'regional-fires',
            img_src: 'aerosols_map_thumbnail.webp',
            alt: 'Map of Alaska showing smoke plume trajectories from regional wildfires to the Juneau Icefield.',
            chartOrder: 2,
            description: 'Smoke plumes from regional wildfires affect glaciers.'
        },
        {
            title: 'Total consumptive use value of recreationally harvested fishes',
            project: 'Rec Fish as Food',
            vizKey: 'FishAsFoodCirclePacking',
            vizRoute: 'rec-fish-consumption-value',
            img_src: 'circle-pack-thumbnail.webp',
            alt: 'Circle packing diagram showing the economic value of recreationally harvested fish, with nested circles representing fish families, species, and countries.',
            chartOrder: 2,
            description: 'The consumption of recreationally harvested inland fish has economic value.'
        },
        {
            title: 'Inland recreational fish consumption',
            project: 'Rec Fish as Food',
            vizKey: 'FishAsFoodSankey',
            vizRoute: 'rec-fish-consumption',
            img_src: 'FishAsFood_sankey_thumbnail.webp',
            alt: 'Sankey diagram showing the flow of inland recreational fish harvest from fish families to species to countries.',
            chartOrder: 1,
            description: 'Recreational harvest of inland fish for consumption is substantial.'
        },
        {
            title: 'Beaufort Sea sediment coring',
            project: 'Beaufort Sea',
            vizKey: 'BeaufortSeaCore',
            vizRoute: 'beaufort-sea-sediment-coring',
            img_src: 'BeaufortSeaCore_thumbnail.webp',
            alt: 'Illustrated sequence showing how USGS scientists collect and process ocean sediment cores from the Beaufort Sea.',
            chartOrder: 1,
            description: 'Sediment cores can help build past and present climates.'
        },  
        {
            title: 'Beaufort Sea species',
            project: 'Beaufort Sea',
            vizKey: 'BeaufortSeaSpecies',
            vizRoute: 'beaufort-sea-species',
            img_src: 'BeaufortSeaSpecies_thumbnail.webp',
            alt: 'Charts showing changes in relative abundance of foraminifera and ostracode microfossil species over the past 2000 years in the Beaufort Sea.',
            chartOrder: 2,
            description: 'Microfossils show that Arctic waters are changing.'
        },   
        {
            title: 'Beaufort Sea timeline',
            project: 'Beaufort Sea',
            vizKey: 'BeaufortSeaTimeline',
            vizRoute: 'beaufort-sea-timeline',
            img_src: 'BeaufortSeaTimeline_thumbnail.webp',
            alt: 'Timeline visualization of climate proxy data from Beaufort Sea sediment cores spanning the past 2000 years.',
            chartOrder: 3,
            description: 'Microfossil records help reconstruct past climates.'
        }   
    ]
};