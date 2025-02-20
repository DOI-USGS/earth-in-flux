export default {
    chartGridItems: [
        //we'll be replacing these img_src with paths to location on s3
        //vizRoutes will direct to appropriate subpage
        {
            title: 'Inland fisheries are threatened',
            project: 'Findex',
            vizKey: 'ThreatSankey',
            vizRoute: 'inland-fish-total-threats',
            img_src: 'findex_sankey_thumbnail.webp',
            alt: '',
            chartOrder: 1,
            description: 'Land use change is threatening inland fisheries.'
        },
        {
            title: 'Explore the Juneau Icefield',
            project: 'Fire in Ice',
            vizKey: 'GlacierScan',
            vizRoute: 'glacier-scan',
            img_src: 'glacial_mri_thumbnail.webp',
            alt: '',
            chartOrder: 1,
            description: 'Glacier ice records change.'
        },
        {
            title: 'Wildfire aerosols',
            project: 'Fire in Ice',
            vizKey: 'WildfireAerosols',
            vizRoute: 'wildfire-aerosols',
            img_src: 'wildfire_aerosols_thumbnail.webp',
            alt: '',
            chartOrder: 2,
            description: 'Wildfire aerosols are preserved in glaciers.'
        },          
        {
            title: 'Regional wildfire deposition',
            project: 'Fire in Ice',
            vizKey: 'RegionalFires',
            vizRoute: 'regional-fires',
            img_src: 'aerosols_map_thumbnail.webp',
            alt: '',
            chartOrder: 2,
            description: 'Smoke plumes from regional wildfires affect glaciers.'
        },
        {
            title: 'Global economic value of recreationally fished species',
            project: 'Fish as Food',
            vizKey: 'FishAsFoodCirclePacking',
            vizRoute: 'inland-rec-fish-value',
            img_src: 'circle-pack-thumbnail.webp',
            alt: '',
            chartOrder: 2,
            description: 'Inland recreational fishing contributes economic value.'
        },
        {
            title: 'Global harvest of recreationally fished species',
            project: 'Fish as Food',
            vizKey: 'FishAsFoodSankey',
            vizRoute: 'inland-rec-fish-harvest',
            img_src: 'FishAsFoodSankey_thumbnail.webp',
            alt: '',
            chartOrder: 1,
            description: 'Inland recreational fishing harvest is substantial.'
        },
        {
            title: 'Beaufort Sea sediment coring',
            project: 'Beaufort Sea',
            vizKey: 'BeaufortSeaCore',
            vizRoute: 'beaufort-sea-sediment-coring',
            img_src: 'BeaufortSeaCore_thumbnail.webp',
            alt: '',
            chartOrder: 1,
            description: 'Sediment cores can help build past and present climates.'
        },  
        {
            title: 'Beaufort Sea species',
            project: 'Beaufort Sea',
            vizKey: 'BeaufortSeaSpecies',
            vizRoute: 'beaufort-sea-species',
            img_src: 'BeaufortSeaSpecies_thumbnail.webp',
            alt: '',
            chartOrder: 2,
            description: 'Microfossils show that Arctic waters are changing.'
        },   
        {
            title: 'Beaufort Sea timeline',
            project: 'Beaufort Sea',
            vizKey: 'BeaufortSeaTimeline',
            vizRoute: 'beaufort-sea-timeline',
            img_src: 'BeaufortSeaTimeline_thumbnail.webp',
            alt: '',
            chartOrder: 3,
            description: 'Microfossil records help reconstruct past climates.'
        }   
    ]
};