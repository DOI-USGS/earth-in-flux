export default {
    chartGridItems: [
        //we'll be replacing these img_src with paths to location on s3
        //vizRoutes will direct to appropriate subpage
        {
            title: 'Inland fisheries are threatened',
            project: 'Findex',
            vizKey: 'ThreatBumpChart',
            vizRoute: 'inland-fish-threats',
            img_src: 'ThreatBumpChart_thumbnail.png',
            alt: '',
            chartOrder: 1,
            description: 'Inland fisheries are threatened.'
        },     
        {
            title: 'Climate vulnerability of recreationally fished inland fish species',
            project: 'Fish as Food',
            vizKey: 'FishAsFoodLinkChart',
            vizRoute: 'inland-rec-fish-climate',
            img_src: 'FishAsFoodLinkChart_thumbnail.png',
            chartOrder: 3,
            alt: '',
            description: 'Inland fish that are recreationally fished are vulnerable'
        },
        {
            title: 'Glacier/Topography 3D Cross-Section Scan',
            project: 'Fire in Ice',
            vizKey: 'GlacierScan',
            vizRoute: 'glacier-scan',
            img_src: 'glacial_mri_thumbnail.png',
            alt: '',
            chartOrder: 1,
            description: 'Ice cores can record changes in wildfire prevalence.'
        },
        {
            title: 'Wildfire Aerosols',
            project: 'Fire in Ice',
            vizKey: 'WildfireAerosols',
            vizRoute: 'wildfire-aerosols',
            img_src: 'placeholder_thumbnail.png',
            alt: '',
            chartOrder: 2,
            description: 'Wildfires are depositing aerosols on glaciers.'
        },          
        {
            title: 'Wildfire Aerosols',
            project: 'Fire in Ice',
            vizKey: 'Aerosols',
            vizRoute: 'aerosol-paths',
            img_src: 'aerosols_thumbnail.png',
            alt: '',
            chartOrder: 3,
            description: 'Wildfire particles are deposited on glaciers.'
        },
        {
            title: 'Global economic value of recreationally fished species',
            project: 'Fish as Food',
            vizKey: 'FishAsFoodCirclePacking',
            vizRoute: 'inland-rec-fish-value',
            img_src: 'circle-pack-thumbnail.png',
            alt: '',
            chartOrder: 2,
            description: 'Inland recreational fishing contributes economic value.'
        },
        {
            title: 'Global harvest of recreationally fished species',
            project: 'Fish as Food',
            vizKey: 'FishAsFoodSankey',
            vizRoute: 'inland-rec-fish-harvest',
            img_src: 'FishAsFoodSankey_thumbnail.png',
            alt: '',
            chartOrder: 1,
            description: 'Inland recreational fishing harvest is substantial.'
        },
        {
            title: 'Beaufort Sea Sediment Coring',
            project: 'Beaufort Sea',
            vizKey: 'BeaufortSeaCore',
            vizRoute: 'beaufort-sea-sediment-coring',
            img_src: 'BeaufortSeaCore_thumbnail.PNG',
            alt: '',
            chartOrder: 1,
            description: 'Sediment cores can help build past and present climates.'
        },  
        {
            title: 'Beaufort Sea Species',
            project: 'Beaufort Sea',
            vizKey: 'BeaufortSeaSpecies',
            vizRoute: 'beaufort-sea-species',
            img_src: 'BeaufortSeaSpecies_thumbnail.png',
            alt: '',
            chartOrder: 2,
            description: 'Microfossils show that Arctic water chemistry is changing.'
        },   
        {
            title: 'Beaufort Sea Timeline',
            project: 'Beaufort Sea',
            vizKey: 'BeaufortSeaTimeline',
            vizRoute: 'beaufort-sea-timeline',
            img_src: 'BeaufortSeaTimeline_thumbnail.PNG',
            alt: '',
            chartOrder: 3,
            description: 'Microfossil records help reconstruct past climates.'
        }   
    ]
};