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
            description: '... and inland fisheries are threatened.'
        },          
        {
            title: 'Glacier/Topography D3 Cross-Section Scan',
            project: 'Fire in Ice',
            vizKey: 'GlacierScan',
            vizRoute: 'glacier-scan',
            img_src: 'glacial_xray_thumbnail.png',
            alt: '',
            description: '...and researchers are studying glacial ice.'
        },       
        {
            title: 'Global economic value of recreationally fished species',
            project: 'Fish as Food',
            vizKey: 'FishAsFoodCirclePacking',
            vizRoute: 'inland-rec-fish-value',
            img_src: 'CirclePacking_thumbnail.png',
            alt: '',
            description: '...and inland recreational fishing contributes economic value.'
        },
        {
            title: 'Beaufort Sea Sediment Coring',
            project: 'Beaufort Sea',
            vizKey: 'BeaufortSeaCore',
            vizRoute: 'beaufort-sea-ice-coring',
            img_src: 'BeaufortSeaCore_thumbnail.PNG',
            alt: '',
            description: '...and sediment cores can help build past and present climates.'
        },  
        {
            title: 'Beaufort Sea Species',
            project: 'Beaufort Sea',
            vizKey: 'BeaufortSeaSpecies',
            vizRoute: 'beaufort-sea-species',
            img_src: 'BeaufortSeaSpecies_thumbnail.png',
            alt: '',
            description: '...and microfossils can be used to indicate these changes.'
        },   
        {
            title: 'Beaufort Sea Timeline',
            project: 'Beaufort Sea',
            vizKey: 'BeaufortSeaTimeline',
            vizRoute: 'beaufort-sea-timeline',
            img_src: 'BeaufortSeaTimeline_thumbnail.PNG',
            alt: '',
            description: '...and communities of microorganisms on the sea floor are affected.'
        }   
    ]
};