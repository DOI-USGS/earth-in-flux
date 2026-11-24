const viteAppAssetURl = import.meta.env.VITE_APP_ASSET_URL;

export default {
    landingPage: {
        pageTitle: "Earth is in flux...",
        pageSubTitle: "USGS scientists research changing environments to inform natural resource management and decision making.",
        collaborationHeading: "The Earth in Flux chart gallery project",
        collaborationDescription: "The Earth in Flux chart gallery is a shared initiative between the U.S. Geological Survey (USGS) <a href='https://www.usgs.gov/mission-areas/water-resources' target='_blank'>Water Resources Mission Area</a> and <a href='https://www.usgs.gov/mission-areas/ecosystems' target='_blank'>Ecosystems Mission Area</a> to communicate key findings of USGS land change science in innovative ways, and to encourage creativity, exploration, and community in data visualization across USGS mission areas.",
        teamHeading: "The USGS Vizlab team",
        teamText: "The USGS Vizlab is a data visualization team within the USGS <a href='https://www.usgs.gov/mission-areas/water-resources' target='_blank'>Water Resources Mission Area</a>. View the Vizlab <a href='https://water.usgs.gov/vizlab/' target='_blank'>portfolio</a>.",
        teamData: [
            { name: "Hayley Corson-Dosch", link: "https://www.usgs.gov/staff-profiles/hayley-corson-dosch", image: "headshots/HCorson-Dosch.png" },
            { name: "Maggie Jaenicke", link: "https://www.usgs.gov/staff-profiles/margaret-maggie-jaenicke", image: "headshots/MaggieJaenicke.png" },
            { name: "Cee Nell", link: "https://www.usgs.gov/staff-profiles/cee-nell", image: "headshots/cee_png.jpeg" },
            { name: "Althea Archer", link: "https://www.usgs.gov/staff-profiles/althea-archer", image: "https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/styles/staff_profile/public/media/images/aaarcher_staff_profile.jpg?h=585bdce6&itok=Z0LQ51Gs" },
            { name: "Jeffrey Kwang", link: "https://www.usgs.gov/staff-profiles/jeffrey-kwang", image: "headshots/jeffrey_kwang_profile.png" },
            { name: "Kaysa Vaarre-Lamoureux", link: "", image: "headshots/kaysa_squarepic.png" },
            { name: "Elmera Azadpour", link: "https://www.usgs.gov/staff-profiles/elmera-azadpour", image: "headshots/elmera_azadpour_square_update.png"}
        ],
        projectsHeading: "USGS land change science projects",
        projectsLeadIn: "The USGS land change science projects currently highlighted in this page are the"
    },
    projects: {
        // keys must match project routes (with '-' replaced with '')
        findex: {
            title: "Findex",
            fullTitle: "Inland Fisheries Threat Index Project",
            logo: "Findex_Logo.png",
            motivation: {
                paragraph1: "Global freshwater biodiversity is declining rapidly. Although rivers, lakes, and other inland waters make up less than 0.01% of Earth's available surface water, they support 51% of all fish species. Inland fisheries also provide food for billions of people and support the livelihoods of millions worldwide.",
                paragraph2: "Despite their importance, inland fishes are some of the most threatened taxa on the planet. Changes to rivers and water flows, damage to streamside habitats, invasive species, climate change, and other pressures have put about one-third of inland fish species at risk of extinction. Yet it is difficult to consistently monitor and asses inland fisheries around the world because they are widely dispersed and often poorly represented in formal markets and data systems.",
                paragraph3: "Here, we present the first global metric of threats to inland fisheries at the river-basin scale. This work combines literature synthesis, expert knowledge, and computational modeling to create a standardized assessment that  serves as a potential risk indicator for freshwater ecosystem status and its capacity to support inland fisheries. These results show that most threats to inland fisheries come from outside the fishery sector, predominately from land use change. This index can help managers identify where risks are greatest and better target limited resources for watershed management, sustainable fisheries, and the people who depend on them. "
            },
            teamText: null,
            teamData: [
                { name: "Gretchen Stokes", link: "", image: "headshots/GretchenStokes.png" },
                { name: "Abigail Lynch", link: "https://www.usgs.gov/staff-profiles/abigail-lynch", image: "https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/styles/staff_profile/public/media/images/abbyprofile.png?h=a35909eb&itok=ob1vliWk" }
            ]
        },
        fireinice: {
            title: "Fire in Ice",
            fullTitle: "Fire in Ice",
            motivation: {
                paragraph1: "Glaciers and ice caps around the world are melting.",
                paragraph2: "Glaciers serve as water towers that store freshwater that is essential for drinking water and agriculture. When glaciers melt, this freshwater is lost to the ocean, leading to sea level rise. Glacial retreat is occurring due to increased temperatures, water collecting under glaciers and hastening their movement, alterations in atmospheric circulation, and the deposition of dark aerosols on ice.",
                paragraph3: "Dark aerosols absorb heat and can increase melting on glacier surfaces. Black carbon, or soot, is a type of dark aerosol that is produced by both industrial burning of fossil fuels and forest fires. Nagorski et al. (2019) determined that black carbon increases the amount of heat absorbed by the Juneau Icefield and accelerates the melting.",
                paragraph4: "The Fire in Ice Project collected a series of ~10-meter snow cores across the Juneau Icefield, AK, which includes some of the most rapidly melting glaciers on the planet and is a major contributor to sea level rise. The goal was to determine if and how wildfire aerosols are affecting the Juneau Icefield by studying material preserved in the core layers.",
                paragraph5: "Scientists analyzed the core for specific sugars (levoglucosan, mannosan, and galactosan) that are only produced by burning vegetation. These sugars are transported in smoke plumes alongside dark aerosols and are also deposited on the icefield. All three sugars were present in quantifiable amounts, demonstrating that wildfires are unequivocally depositing aerosols onto the icefield. Ratios between the sugars can provide information on the vegetation type—hardwoods or softwoods—that burned in the fires. The presence of hardwood fires, coupled with satellite data of smoke transport, suggests that fires as far away as East Asia may affect the Juneau Icefield."
            },
            teamText: null,
            teamData: [
                { name: "Natalie Kehrwald", link: "https://www.usgs.gov/staff-profiles/natalie-m-kehrwald", image: "https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/styles/staff_profile/public/media/images/NatalieKehrwald_cropped.jpg?itok=oLP9Dl1H" },
                { name: "Jeramy R. Jasmann", link: "https://www.usgs.gov/staff-profiles/jeramy-r-jasmann", image: "https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/styles/staff_profile/public/thumbnails/image/Jasmann_Jeramy_USGS.JPG?itok=ADC3fg6X" },
                { name: "Michelle Leung", link: "https://www.usgs.gov/staff-profiles/michelle-leung", image: "https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/styles/staff_profile/public/media/images/MichelleLeung.JPG?h=bcaad0e8&itok=LXDUuI9W" },
                { name: "Morgan Schachterle", link: "https://www.usgs.gov/staff-profiles/morgan-schachterle", image: "https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/styles/staff_profile/public/media/images/Morgan%20_Schachterle_profileImage.jpeg?itok=r244L8u8" }
            ]
        },
        recfishasfood: {
            title: "Rec Fish as Food",
            fullTitle: "Recreational Inland Fisheries as Food",
            pageHasSubtitle: true,
            pageSubtitle: "Inland recreational fisheries contribute nutritional benefits and economic value but are vulnerable to climate change",
            motivation: {
                paragraph1: "Inland recreational fishing, defined as primarily leisure-driven fishing in freshwaters, is a popular past-time that can provide substantial contributions to human consumption. However, these benefits to food systems and food security are often overlooked at global scales. Here, we established a baseline of inland recreational consumption estimates by fish species and country to identify the nutritional composition, consumptive use value, and climate vulnerability of this recreational consumption.",
                paragraph2: "To learn more, read the USGS Science Spotlight <a href='https://www.usgs.gov/programs/climate-adaptation-science-centers/news/recreational-inland-fishing-fun-and-food' target='_blank'>blog post</a> highlighting this research, visit the <a href='https://cascprojects.org/#/project/5050cb0ee4b0be20bb30eac0/6026e1a9d34eb120311394c4' target='_blank'>project page</a>, or explore the data on the <a href='https://rconnect.usgs.gov/recfishfood/' target='_blank'>Recreational Inland Fisheries as Food web app</a>."
            },
            teamText: null,
            teamData: [
                { name: "Holly Embke", link: "https://www.usgs.gov/staff-profiles/holly-embke", image: "https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/styles/staff_profile/public/media/images/embke.jpg?h=f1072bca&itok=xyG6frc8" },
                { name: "Abigail Lynch", link: "https://www.usgs.gov/staff-profiles/abigail-lynch", image: "https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/styles/staff_profile/public/media/images/abbyprofile.png?h=a35909eb&itok=ob1vliWk" }
            ]
            
        },
        beaufortsea: {
            title: "Beaufort Sea",
            fullTitle: "Beaufort Sea",
            motivation: {
                paragraph1: "The Arctic Ocean is undergoing dramatic sea ice reduction and warming conditions.",
                paragraph2: "These changes affect the sealife of the region, including bottom-dwelling organisms and the marine mammals, seabirds, and fish that rely on them for food. The researchers of this project use microfossils from sediment cores taken in the Beaufort Sea to rebuild the climate patterns, sea ice and circulation, and ecosystems from the past 2000 years.",
                paragraph3: "By studying bottom-dwelling microorganisms like microscopic arthropods (called Ostracodes) and single-celled protists (called Foraminifera), the Beaufort Sea researchers can compare the historical conditions of the Arctic and better understand the effects of recent climate change in the region. Read more about this USGS science <a href='https://www.usgs.gov/programs/ecosystems-land-change-science-program/science/land-sea-linkages-arctic' target='_blank'>here</a>."
            },
            teamText: "Laura Gemery is an ecologist, and Julia Seidenstein, Jason Addison and Thomas Cronin are geologists within the <a href='https://www.usgs.gov/programs/ecosystems-land-change-science-program' target='_blank'>Ecosystems Land Change Science program</a> of the USGS Ecosystems Mission area. Each specialize in analyzing proxies in sediment cores.",
            teamData: [
                { name: "Laura Gemery", link: "https://www.usgs.gov/staff-profiles/laura-gemery", image: "headshots/l_Gemery mustang Oden.png" },
                { name: "Julia Seidenstein", link: "https://www.usgs.gov/staff-profiles/julia-seidenstein", image: "headshots/Julia Seidenstein_professional pic.png" },
                { name: "Jason Addison", link: "https://www.usgs.gov/staff-profiles/jason-addison", image: "headshots/Jason Addison.png" },
                { name: "Thomas Cronin", link: "https://www.usgs.gov/staff-profiles/thomas-cronin", image: "headshots/T_Cronin head shot.png" }
            ]
        }
    },
    visualizations: {
        GlacierScan: {
            paragraph1: "The Juneau Icefield is just north of Juneau, Alaska and extends into Canada. In 2016 and 2017, researchers collected snow cores across the icefield to determine if the ice traps records of wildfires.",
            paragraph1Mobile: "The Juneau Icefield is just north of Juneau, Alaska. In 2016 and 2017, researchers collected snow cores across the icefield to determine if the ice traps records of wildfires.",
            promptDesktop: "Hover over the map to explore the topography and learn about coring a glacier.",
            promptMobile: "Tap the map to explore the topography and learn about coring a glacier.",
            heading: "How are ice cores collected?",
            paragraph2: "The ideal location for drilling a snow core to identify biomarkers from wildfire events is the highest, flattest section of a glacier, as this region minimizes any influence from the glacier flow. Snow builds up layer by layer on a glacier surface and eventually compresses into ice. High, flat drilling locations increase the possibility that packed layers remain horizontal.",
            paragraph3: "Once researchers have arrived at the drilling location, the first thing that they do is create some sort of shelter for comfort in case the weather becomes bad. For short (~10 m) snow cores that can be drilled in a day, this shelter is often a snow pit with a bench and a tarp to block wind. Once the pit is dug, scientists can start drilling.",
            paragraph4: "Snow core drills either have a cutting head (like a Felix corer) or cutting teeth (like a Kovacs corer). The teeth cut the packed snow. Parts called ‘core dogs’ grab the core so that you can retrieve it from the bottom of the drilled hole. The drill itself is a metal cylinder and has threads on the outside that help excavate the cut core, allowing the drill to keep cutting deeper.",
            paragraph5: "Researchers attach a handle to the top and turn the drill until they have collected a meter of packed snow. This meter of snow is then passed to the people in the snow pit who are processing the core. This processing entails documenting the snow stratigraphy, measuring the core, and weighing the section to determine its mass.",
            paragraph6: "The drilling and processing continue simultaneously, with team members measuring one core section while the next section is being drilled. The snow core drill returns to the initial hole and retrieves the next meter of packed snow. As the drill goes deeper and deeper into the glacier, metal extension rods are connected to the top of the drill to extend the drill length until eventually reaching ~10 meters below the surface.",
            corerAlt1: "Top view of a Kovacs snow core drill. It is a long narrow cylinder, with cutting teeth at one end and threads that run up the length ofo the cylinder to a black t-shaped handle.",
            corerAlt2: "Close-up view of the core dogs on a Kovacs snow core drill, which grab the core so that it can be retrieved.",
            corerAlt3: "Close-up view of the cutting head on a Kovacs snow core drill.",
            photo010: "The Juneau Icefield Research Program (JIRP) established a series of camps across the icefield that allow students and faculty the opportunity to live and work on the icefield throughout the summer program. Clear, sunny days on the Juneau Icefield are rare and much appreciated, and researchers spend every possible moment outside enjoying the spectacular scenery.",
            photo010Mobile: "Clear, sunny days on the Juneau Icefield are rare and much appreciated, and researchers spend every possible moment outside enjoying the scenery.",
            photo010Alt: "A small building is perched on an outcrop of rocks above the vast expanse of the Juneau Icefield. It is a sunny day with clear blue skies, and scientists stand outside the building taking in the view.",
            photo018: "Drilling and processing snow cores involves teamwork. While a new core section is being drilled, team members measure physical properties (e.g., volume, mass) and stratigraphy (i.e., the ice layering) of the previous core section. Processing the core in a snow pit allows shelter from the wind and creates a stable bench for examining the cores.",
            photo018Mobile: "While a new core section is being drilled, team members measure physical properties (e.g., volume, mass) and stratigraphy (i.e., the ice layering) of the previous core section.",
            photo018Alt: "Five scientists work to collect snow cores at a site in the middle of the vast icefield. One scientist is standing and turning a snow corer to collect a sample. Two scientists in a adjacent pit work to process the core as two scientists sit on the ground taking notes.",
            photo051: "Scientist sample seasonal snow and firn (snow older than one year) on Lemon Creek Glacier, measuring density, noting any layers of interest, and collecting samples for chemical analysis. The coring device can obtain one meter of snow and firn at a time, and by adding extensions to the drill, it is possible to carefully sample many meters deep without having to dig an enormous hole. Team members take careful notes to ensure records are complete.",
            photo051Mobile: "The coring device can obtain one meter of snow and firn at a time, and by adding extensions to the drill, it is possible to carefully sample many meters deep without having to dig an enormous hole.",
            photo051Alt: "Four scientists stand on a plain of snow atop the Lemon Creek Glacier. One scientist is holding a tall vertical rod, at the base of which is a snow corer that the scientists are using to collect samples of the snow.",
            photo085: "Sunsets on the Juneau Icefield are spectacular, with the combination of jagged peaks and low light gleaming on the snow.",
            photo085Mobile: "Sunsets on the Juneau Icefield are spectacular, with the combination of jagged peaks and low light gleaming on the snow.",
            photo085Alt: "A view of the sun setting behind a jagged range of mountain peaks. The light from the sun is reflecting off of a snowfield in the foreground.",
            photo138: "Travel between sites on the icefield involves a combination of skiing across glaciers and hiking over slopes. In order to be able to traverse to a new camp, the teams must be well-practiced in safety skills such as crevasse rescue and rope management.",
            photo138Mobile: "Travel between sites on the icefield involves a combination of skiing across glaciers and hiking over slopes. Teams use safety skills like crevasse rescue and rope management.",
            photo138Alt: "Three scientists traversing the icefield carry large packs that have skis strapped onto them.",
            photo140: "Natalie Kehrwald takes a break from drilling snow cores on the Taku Glacier, one of the deepest temperate mountain glaciers in the world.",
            photo140Mobile: "Natalie Kehrwald takes a break from drilling snow cores on the Taku Glacier, one of the deepest temperate mountain glaciers in the world.",
            photo140Alt: "USGS scientist Natalie Kehrwald smiles as she stands facing the camera, wearing a magenta jacket, a black hat, dark sunglasses, and a safety harness. The icefield extends into the distance behind her.",
            photo156: "Sunny summer days create a sun-cupped surface texture on the snow surface. The passage of scientists on skis and snowmobiles is visible here, showing a trail from camp to the coring site.",
            photo156Mobile: "Sunny summer days create a sun-cupped surface texture on the snow surface. The passage of scientists on skis and snowmobiles is visible here.",
            photo156Alt: "Looking out across the icefield, the surface is of the snow appears dimpled. A strip of snow that extends into the distance has been visibly smoothed and packed down by scientists skiing and driving snowmobiles.",
            photo203: "Glaciers scour the landscape, creating U-shaped valleys between surrounding peaks. At their terminus (where they end), glacier melt creates lakes and rivers.",
            photo203Mobile: "Glaciers scour the landscape, creating U-shaped valleys between surrounding peaks. Where they end, glacier melt creates lakes and rivers.",
            photo203Alt: "A view of the jagged end of a glacier where it terminates into a glacial lake.",
            photo021: "The movement of glaciers pushes sediments both in front of and to the sides of the ice. The mounds of sediments are called moraines and remain even after the ice has melted. Moraines provide evidence of past glacier activity in locations that are now free of ice.",
            photo021Mobile: "The movement of glaciers pushes sediments both in front of and to the sides of the ice. The mounds of sediments are called moraines.",
            photo021Alt: "From the high vantage point of a ridge, you can see the glacier snaking down a deep u-shaped valley that is surrounded by jagged peaks. On the surface of the glacier, running parallel to the sides of the valley, are two sinous dark bands of material that has collected on the surface of the glacier as it has moved."
        },
        RegionalFires: {
            paragraph1: `From analysis of <a href='/${viteAppAssetURl}/earth-in-flux/fire-in-ice/glacier-scan' target='_blank'>glacier snow cores</a>, we know that wildfires burning softwoods like pine, fir, and spruce have <a href='/${viteAppAssetURl}/earth-in-flux/fire-in-ice/wildfire-aerosols' target='_blank'>deposited aerosols</a> on the Juneau Icefield. We also know that Alaska has many softwood forests, and that some have burned. But how can we tell which regional fires deposited aerosols on the glacier?`,
            paragraph2: "Researchers use an atmospheric model (<a href='https://www.ready.noaa.gov/HYSPLIT.php' target='_blank'>HYSPLIT</a>) to trace the potential path of smoke particles generated by known wildfires, identifying 'candidate' fires that could be the source of aerosols deposited on the icefield during the sampled period."
        },
        BeaufortSeaCore: {
            heading1: "Why collect ocean sediment cores?",
            intro1: "Ocean sediments are one of the best archives of past ocean conditions and changes to the climate throughout Earth's history.",
            intro2: "USGS scientists collect sediment cores from the ocean floor to reconstruct past environmental conditions. This knowledge helps us to understand processes that influence natural climate variability.",
            heading2: "How do USGS scientists collect ocean sediment core records?",
            paragraph1: `Collecting sediment core records from the ocean floor begins with USGS scientists traveling to their study site in Utqiagvik, Alaska. At this site, scientists used recent climatological measurements and several types of climate ‘proxies’ to reconstruct the climate history of the Beaufort Sea and to better understand recent climate change in the Arctic Ocean. Climate proxies in ocean sediment cores can include preserved physical or chemical properties or biological organisms, like algae or <a href='/${viteAppAssetURl}/earth-in-flux/beaufort-sea/beaufort-sea-species' target='_blank'>microorganisms with shells</a>. These proxies serve as indicators of ocean- and climate-related conditions at the time the sediment was deposited.`,
            paragraph2: "Once they arrived in Alaska, the researchers boarded the U.S. Coast Guard icebreaker <i>Healy</i> and traveled to sampling locations along the continental shelf, where the Mackenzie River flows into the Beaufort Sea.",
            paragraph3: "Researchers collected sediment cores by deploying a weighted piston-driven coring device from a platform on the stern of the ship. The device retrieved a vertical core of soft sediments from the bottom of the ocean. Because sediment accumulates on the seafloor over time, the sediment at the top of the core is younger than sediment at the bottom.",
            paragraph4: `Onboard the ship, researchers cut the sediment core into one- and five-centimeter-thick slices for analysis, representing ~4- and ~24-year periods, respectively. The total time span that this particular sediment core represents is about <a href='/${viteAppAssetURl}/earth-in-flux/beaufort-sea/beaufort-sea-timeline' target='_blank'>2000 years</a>. The timeline of sediment deposition was established by measuring <span class='tooltip-group'><span class='tooltip-span'>Cesium (Cs)</span><span id='cesium-tooltip' class='tooltiptext'>Cesium-137 dating measures a radioactive isotope released during nuclear tests in the 1950s and 1960s, which indicates when sediment layers were deposited from the mid-20th century onward.</span></span> and <span class='tooltip-group'><span class='tooltip-span'>Lead (Pb)</span><span id='lead-tooltip' class='tooltiptext'>Lead-210, a decay product of radon originating from uranium in the soil, accumulates in sediments over time. Measuring its concentration helps estimate the age of sediment layers over the last ~150 years.</span></span> radioisotopes and by <span class='tooltip-group'><span class='tooltip-span'>radiocarbon dating</span><span id='radiocarbon-tooltip' class='tooltiptext'>Researchers use radiocarbon dating to extend the record further back in time (up to ~50,000 years). This method relies on carbon-14, which is produced in the atmosphere and absorbed by living organisms. Once an organism dies, it stops taking in Carbon-14, which then decays at a known rate.</span></span> fossilized shells down the length of the core.`,
            paragraph5: "After processing the core, the researchers identified all the microfossils in each slice. The microfossils ranged in size from 0.5 to 2 millimeters. Most of them were a little larger than the period at the end of this sentence.",
            paragraph6: `Scientists used a damp brush to move the microfossils from the sediment and place them on a slide for identification under a microscope. The microfossils examined in this study include <a href='/${viteAppAssetURl}/earth-in-flux/beaufort-sea/beaufort-sea-species' target='_blank'>foraminifera and ostracoda</a>, which are excellent proxy indicators for past environments because their presence and growth were influenced by specific climate and ocean conditions when they were alive.",
            alt1: "Illustrated world map showing a cartoon airplane flying from Reston, Virginia at USGS headquarters to the Beaufort Sea at the northeast corner of Alaska.`,
            alt2: "Zoomed in illustrated map of Alaska and northwest Canada. A cartoon icebreaker ship is offshore from the mouth of the Mackenzie River where it opens up into the Beaufort Sea. Text labels Alaska, Canada, Mackenzie River, and Beaufort Sea.",
            alt3: "Animated cartoon gif showing the icebreaker ship pushing through some ice to get to the study site. The landscape shows ice-covered Beaufort Sea with a dull, cold blue-grey sky and distant clouds.",
            alt4: "Two side-by-side photographs of the ship on the left and people on deck of the ship on the right. The images are composed in a scrapbook-style alongside a cartoon beluga whale. The ship is large and bright red with the U.S. Coast Guard label on the side. Behind the ship are mountains covered in ice. The people on the deck of the ship are watching the sediment core as it enters the ice. They're wearing bright red Coast Guard jackets to stay safe and warm.",
            alt5: "Animated cartoon gif showing the Coast Guard ship on top of the icy sea and dropping the sediment core down into the ocean. When the core gets to the bottom, it removes a core of the ocean floor and then is returned to the ship with the sediment sample. Swimming under the sea is a cartoon beluga whale.",
            alt6: "Three photographs of the researchers collecting the core, processing it, and then also showing the core once it is in the lab. The images are composed in a scrapbook-style alongside the cartoon Coast Guard ship.",
            alt7: "Animated cartoon gif showing the sediment core after it's been pulled up from the ocean floor. As the animation moves forward, the core is sliced into small discs. One example disc is shown zoomed-in with lots of little specks of color, representing the microfossils. The five focal species are also shown zoomed-in near the sample to represent that they are identified from these sediment discs.",
            alt8: "Three photographs of the researchers examining and identifying the microfossils. One image shows Laura Gemery, USGS researcher, studying microfossils through a microscope. The middle image shows a zoomed in image of a petri dish and Laura holding a small paint brush and one microfossil, which is barely visible next to the paint brush bristles. The third image shows a close-up of the microscope with some samples laid out underneath it. The scrapbook-style images are surrounded by actual scanning electron images of the different microfossil species."
        },
        BeaufortSeaSpecies: {
            heading1: "The Arctic Ocean is changing",
            paragraph1: "In the last 25 to 50 years, major shifts have taken place in ecological communities in the Beaufort Sea, including within communities of microorganisms. Because of their sensitivity to water chemistry, microorganisms are proxies for local climatic conditions. Taken together with other climate proxies analyzed at this site, the shifts in microorganism abundance indicate that the climate and ocean conditions have recently changed in this region of the Beaufort Sea.",
            heading2: "What can microorganisms tell us about change in the Arctic Ocean?",
            paragraph2: `Ocean floor <a href='/${viteAppAssetURl}/earth-in-flux/beaufort-sea/beaufort-sea-sediment-coring' target='_blank'>sediment cores</a> provide a record of changes in microorganism communities. Microorganisms like foraminifera and ostracodes are identified in the sediment by their shells, which can take many different forms. These organisms build their shells using calcium carbonate and other elements taken directly from ocean water. Under the right conditions, these shells can become fossilized in sediments, creating microfossils and preserving a long-term record of past ocean conditions. Ocean water pH, salinity, turbidity, and temperature affect the sturdiness of the shells and the availability of different food sources. As a result, the presence and abundance of individual species is closely tied to environmental conditions. Because of the known ecological preferences of each species, microfossils act as proxies of past ocean conditions and as indicators of environmental change in the Arctic Ocean.`,
            subheading1: "Foraminifera",
            paragraph3: "Foraminifera are single-celled protists. Foraminifera secrete a shell called a 'test' made of calcium carbonate. Because each species creates tests with unique structures and shapes, researchers can identify individual species in the sedimentary fossil record.",
            subheading2: "Ostracodes",
            paragraph4: "Ostracodes are crustaceans that secrete calcium carbonate to form a bivalve shell (like that of a clam).",
            tabData: [
                {
                    tabTitle: "Spiroplectammina biformis",
                    tabContentTitle: "Spiroplectammina biformis",
                    tabContentTitleID: "spiroplectimmina",
                    tabSpeciesClass: "Foraminifera",
                    tabText: "The relative abundance of <span class='scientificName'>Spiroplectammina biformis</span> has been higher in the past few decades than in the last 2000 years. Species within the <span class='scientificName'>Spiroplectammina</span> genus create tests that are 'agglutinated,' or formed by cementing together particles of sand and sediment. Agglutinated species are able to withstand harsher environments than other types of foraminifera. High abundances of agglutinated species often indicate corrosive (acidic) and turbid (cloudy) conditions. The recent increase in the relative abundance of <span class='scientificName'>Spiroplectammina biformis</span> suggests that eroded permafrost carried into the Beaufort Sea by the Mackenzie River as a consequence of climate change may be making bottom ocean waters more acidic and less hospitable to other species that have carbonate shells.",
                    tabContentImageSuffix: "2c",
                    tabPrefixImageName: "F_Spiroplectammina",
                    tabImageAlt: "Scatterplot of the relative abundance of Spiroplectammina biformis over the past 2000 years. The relative abundance was low until about 1000 C.E., then it had some periods of being a bit higher (about 25% to 40%) until about 1500 C.E.. In the last 100 years or so, the relative abundance has been increasing rapidly and is above 50% relative abundance. A cartoon image of the microfossil is overlaid in the upper left side of the plot."
                },
                {
                    tabTitle: "Cassidulina reniforme",
                    tabContentTitle: "Cassidulina reniforme",
                    tabContentTitleID: "cassidulina",
                    tabSpeciesClass: "Foraminifera",
                    tabText: "The relative abundance of <span class='scientificName'>Cassidulina reniforme</span> has declined very recently. This shift may be due to the increased abundance of agglutinated species like <span class='scientificName'>Spiroplectammina biformis</span>, which can survive in more corrosive conditions.",
                    tabContentImageSuffix: "2a",
                    tabPrefixImageName: "F_Cassidulina",
                    tabImageAlt: "Scatterplot of the relative abundance of Cassidulina reniforme over the past 2000 years. The relative abundance was fairly stable, oscillating around about 40% relative abundance until the last 500 years or so when the abundance has declined to a low recently near 0% abundance. A cartoon image of the microfossil is overlaid in the upper left side of the plot."
                },
                {
                    tabTitle: "Elphidium excavatum",
                    tabContentTitle: "Elphidium excavatum",
                    tabContentTitleID: "elphidium",
                    tabSpeciesClass: "Foraminifera",
                    tabText: "Like <span class='scientificName'>Cassidulina reniforme</span>, <span class='scientificName'>Elphidium excavatum</span> has declined in relative abundance in the Beaufort Sea in recent years. As climate change (and potentially an increased sediment load in the Mackenzie River) makes local ocean water more acidic, the aquatic environment becomes less favorable to organisms that build calcium carbonate shells.",
                    tabContentImageSuffix: "2b",
                    tabPrefixImageName: "F_Elphidium",
                    tabImageAlt: "Scatterplot of the relative abundance of Elphidium excavatum over the past 2000 years. The relative abundance was fairly stable, oscillating around about 25% relative abundance until the last 500 years or so when the abundance has declined to a low recently near 0% abundance. A cartoon image of the microfossil is overlaid in the upper left side of the plot."
                },
                {
                    tabTitle: "Kotoracythere arctoborealis",
                    tabContentTitle: "Kotoracythere arctoborealis",
                    tabContentTitleID: "kotorachythere",
                    tabSpeciesClass: "Ostracode",
                    tabText: "<span class='scientificName'>Kotoracythere arctoborealis</span> is a species of ostracode that is not as common and is less widespread than other species of ostracode in the Arctic Ocean. Its ecological preferences are more limiting, and it prefers protected areas, like bays with comparatively warmer temperatures. In the Beaufort Sea, <span class='scientificName'>Kotoracythere arctoborealis</span> has had relatively low, but steady, abundance until recently, when it increased in relative abundance. This shift is likely due to changes in ocean water salinity and warmer summer water temperatures that favor <span class='scientificName'>Kotoracythere arctoborealis</span>.",
                    tabContentImageSuffix: "3a",
                    tabPrefixImageName: "O_Kotoracythere",
                    tabImageAlt: "Scatterplot of the relative abundance of Kotoracythere arctoborealis over the past 2000 years. The relative abundance was fairly stable and low (about 10% or less) over much of the timespan. In the last 500 years or so, the abundance has increased to almost 25%. A cartoon image of the microfossil is overlaid in the upper left side of the plot."
                },
                {
                    tabTitle: "Paracyprideis pseudopunctillata",
                    tabContentTitle: "Paracyprideis pseudopunctillata",
                    tabContentTitleID: "paracyprideis",
                    tabSpeciesClass: "Ostracode",
                    tabText: "<span class='scientificName'>Paracyprideis pseudopunctillata</span> is common in polar regions with bottom waters that are very cold year-round, and historically has been one of the most dominant species in the Beaufort Sea. It is becoming less dominant as conditions favor other species like <span class='scientificName'>Kotocythere arctoborealis</span> and <span class='scientificName'>Spiroplectammina biformis</span>.",
                    tabContentImageSuffix: "3b",
                    tabPrefixImageName: "O_Paracyprideis",
                    tabImageAlt: "Scatterplot of the relative abundance of Paracyprideis pseudopunctillata over the past 2000 years. The relative abundance was fairly stable, oscillating around about 60% relative abundance until around 1500 C.E. when it declined to nearly 0%. Then, it increased again back to about 60% abundance until about 1750 C.E., and it has been declining since. The relative abundance was just above 25% for this species in the last sample in 2000 C.E. A cartoon image of the microfossil is overlaid in the upper left side of the plot."
                }
            ],
            heading3: "How are these data collected?",
            paragraph5: "These data are from the mouth of the Mackenzie River where it spills into the Beaufort Sea north of Yukon, Canada. In this dynamic Arctic environment, sea ice extent, ocean water temperature and salinity, and the availability of food can change spatially and temporally. The species living in bottom ocean waters shift in response to these variable conditions. These changes are recorded in the sediment when the microorganisms are preserved as microfossils.",
            paragraph6: `USGS researchers collect <a href='/${viteAppAssetURl}/earth-in-flux/beaufort-sea/beaufort-sea-sediment-coring' target='_blank'>sediment cores</a> and analyze the relative abundance of species within each sediment layer. See the full 2000-year timeline of microfossil composition on the <a href='/${viteAppAssetURl}/earth-in-flux/beaufort-sea/beaufort-sea-timeline' target='_blank'>Beaufort Sea timeline</a> page.`,
        },
        BeaufortSeaTimeline: {
            heading1: "Microfossils help reconstruct climate records",
            paragraph1: `Although climate change is a global phenomenon, it has uneven local consequences. Research into the history of the Earth’s climate provides a deeper understanding of processes that influence natural climate variability. Reconstructed climate records, built using indicators like <a href='/${viteAppAssetURl}/earth-in-flux/beaufort-sea/beaufort-sea-species' target='_blank'>microfossils</a> in ocean floor <a href='/${viteAppAssetURl}/earth-in-flux/beaufort-sea/beaufort-sea-sediment-coring' target='_blank'>sediment cores</a>, can extend modern information back in time. These records provide baselines for natural conditions before humans began modifying Earth processes. In the Arctic, records of past environmental conditions can place current anthropogenic climate warming and sea-ice loss in a long-term context.`,
            paragraph2: "The sediment record captures microfossil species assemblages in the Beaufort Sea over the past 2000 years. The size of each bubble in the chart below is scaled to represent the relative abundance of an individual species of microfossil, including indicator species in the genera <span class='highlight scientificName' id='elphidium'>Elphidium</span>, <span class='highlight scientificName' id='cassidulina'>Cassidulina</span>, <span class='highlight scientificName' id='paracyprideis'>Paracyprideis</span>, <span class='highlight  scientificName' id='kotorachythere'>Kotoracythere</span> and <span class='highlight scientificName' id='spiroplectimmina'>Spiroplectimmina</span> and <span class='highlight scientificName' id='other-species'>other species</span> of ostracodes and foraminifera. The bubble chart shows species assemblages over the 2000-year record. Hover over the bar chart to look at relative abundances of the species within each 100-year window.",
            heading2: "What does this record tell us?",
            paragraph3: "This sedimentary record recounts the history of climate variability and change in the Arctic Ocean over the past 2000 years. Shifts in microorganism abundance, alongside chemical markers like biogenic silica and ratios of carbon and oxygen isotopes (not shown), reflect variability in the salinity and temperature of ocean water, ecological productivity, and the delivery of terrestrial sediment from the Mackenzie River. The record shows summer warming oscillations during the Medieval Climate Anomaly (beginning around 950 CE), and a colder more variable climate from ~1250 to 1900 CE. In the last 60 years, the record reflects larger-scale changes in the Beaufort Sea. Loss of sea ice has led to longer periods of open water. Permafrost melt has increased inflow from the Mackenzie River and made that inflow more acidic. Near the sea floor, conditions are cloudier, and productivity has increased.",
            paragraph4: "Even though shallow continental shelves make up half of the Arctic Ocean, sediment records from high-latitude continental shelves are rare. With this record of climate variability and ecosystem change, scientists can better understand and model climate change in coastal regions of the Arctic.",
            heading3: "How was this record reconstructed?",
            paragraph5: `USGS scientists reconstructed this record using <a href='/${viteAppAssetURl}/earth-in-flux/beaufort-sea/beaufort-sea-sediment-coring' target='_blank'>sediment cores</a> collected from the ocean floor on the Beaufort Sea continental shelf, north of Yukon, Canada. After cutting the sediment cores into slices representing ~4 to 24 years of time, researchers used a microscope to identify the <a href='/${viteAppAssetURl}/earth-in-flux/beaufort-sea/beaufort-sea-species' target='_blank'>species of microfossils</a> in each slice.`,
        },
        FishAsFoodCirclePacking: {
            paragraph1: `Explore the global <span class='tooltip-group'><span class='tooltip-span'>total consumptive use value</span><span id='tcuv-tooltip' class='tooltiptext'>Calculated by multiplying the total kilograms of species-specific edible portions by the shadow prices for each species per kilogram, in U.S. dollars. Species- and country-specific price data were collected from November 2021 to February 2022.</span></span> in U.S. dollars (2022) of <a href='/${viteAppAssetURl}/earth-in-flux/rec-fish-as-food/rec-fish-consumption' target='_blank'>recreationally harvested inland fish species</a>. Click on the largest circles in the diagram to see the consumptive use value of species in each fish family. Click on the nested blue-gray circles to see the value for a species in each country where it is recreationally harvested for consumption. Note that species with consumptive use values under $500,000 are grouped into 'Other' unless only one exists. Some taxonomic data are not available to the species level. Click anywhere outside the circle to zoom out.`,
            familyInfo: {
                "Cyprinidae": {
                    image: "images/common-carp.jpeg",
                    caption: "Common carp <em>(Cyprinus carpio)</em>"
                },
                "Salmonidae": {
                    image: "images/pink-salmon.jpeg",
                    caption: "Pink salmon <em>(Oncorhynchus gorbuscha)</em>"
                },
                "Centrarchidae": {
                    image: "images/black-crappie.jpeg",
                    caption: "Black crappie <em>(Pomoxis nigromaculatus)</em>"
                },
                "Anguillidae": {
                    image: "images/european-eel.jpeg",
                    caption: "European eel <em>(Anguilla anguilla)</em>"
                },
                "Anostomidae": {
                    image: "images/boga.jpeg",
                    caption: "Boga <em>(Megaleporinus obtusidens)</em>"
                },
                "Atherinopsidae": {
                    image: "images/silverside.jpeg",
                    caption: "Argentinian silverside <em>(Odontesthes bonariensis)</em>"
                },
                "Channidae": {
                    image: "images/snakehead.jpeg",
                    caption: "Northern snakehead <em>(Channa argus)</em>"
                },
                "Characidae": {
                    image: "images/hemibrycon.jpeg",
                    caption: "Hemibrycon <em>(Hemibrycon boquiae)</em>"
                },
                "Cichlidae": {
                    image: "images/ zebra-mbuna.jpeg",
                    caption: "Zebra mbuna <em>(Maylandia zebra)</em>"
                },
                "Esocidae": {
                    image: "images/northern-pike.jpeg",
                    caption: "Northern pike <em>(Esox lucius)</em>"
                },
                "Percidae": {
                    image: "images/walleye.jpeg",
                    caption: "Walleye <em>(Sander vitreus)</em>"
                },
                "Siluridae": {
                    image: "images/european-catfish.jpeg",
                    caption: "European catfish <em>(Silurus glanis)</em>"
                },
                "Ictaluridae": {
                    image: "images/catfish.jpeg",
                    caption: "Catfish <em>(Ictalurus punctatus)</em>"
                },
                "Latidae": {
                    image: "images/nile-perch.jpeg",
                    caption: "Nile perch <em>(Lates niloticus)</em>"
                },
                "Erythrinidae": {
                    image: "images/trahira.jpeg",
                    caption: "Trahira <em>(Hoplias malabaricus)</em>"
                },
                "Lotidae": {
                    image: "images/burbot.jpeg",
                    caption: "Burbot <em>(Lota lota)</em>"
                },
                "Osmeridae": {
                    image: "images/european-smelt.jpeg",
                    caption: "European smelt <em>(Osmerus eperlanus)</em>"
                },
                "Percichthyidae": {
                    image: "images/Percichthyidae.jpeg",
                    caption: "Macquarie perch <em>(Macquaria australasica)</em>"
                },
                "Bagridae": {
                    image: "images/Bagridae.jpeg",
                    caption: "Striped dwarf catfish <em>(Mystus vittatus)</em>"
                },
                "Centropomidae": {
                    image: "images/Centropomidae.jpeg",
                    caption: "Common snook <em>(Centropomus undecimalis)</em>"
                },
                "Osteichthyes": {
                    image: "images/squalius_cephalus.png",
                    caption: "European chub <em>(Squalius cephalus)</em>"
                },
                "Pimelodidae": {
                    image: "images/redtail-catfish.jpeg",
                    caption: "Redtail catfish <em>(Phractocephalus hemioliopterus)</em>"
                },
                "Plecoglossidae": {
                    image: "images/Plecoglossidae.jpeg",
                    caption: "Ayu sweetfish <em>(Plecoglossus altivelis)</em>"
                },
                "Serrasalmidae": {
                    image: "images/silver-dollar.jpeg",
                    caption: "Silver dollar <em>(Metynnis luna)</em>"
                },
                "Sinipercidae": {
                    image: "images/mandarin-fish-edit.jpeg",
                    caption: "Mandarin fish <em>(Siniperca chuatsi)</em>"
                }

            },
            defaultFamily: {
                image: "images/default-fish.jpeg",
                text: "Click the circles in the chart to learn more about the consumptive use value of each fish family and species.",
                caption: "Brown trout <em>(Salmo trutta)</em>"
            },
            accordionData: [
                {
                    heading: "Attributions for fish family silhouette images",
                    content: [
                        {
                            type: "text",
                            content: "The silhouette images used for each fish family are available under Creative Commons licenses on <a href='https://www.phylopic.org/' target='_blank'>PhyloPic</a>:<span><ul><li>The default silhouette shown on page load is <a href='https://www.phylopic.org/images/23a7d09d-4a4d-4ad5-ad07-49a6b59a7fba/salmo-trutta' target='_blank'><em>Salmo trutta</em></a> by Carlos Cano-Barbacil, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Anguillidae: <a href='https://www.phylopic.org/images/9ce911e2-a8fa-47a2-afe9-ca18d71be30f/anguilla-anguilla' target='_blank'><em>Anguilla anguilla</em></a> by Ingo Braasch, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Anostomidae: <a href='https://www.phylopic.org/images/73b3e3c5-198a-48c7-8fcf-1ae704c416af/megaleporinus-muyscorum' target='_blank'><em>Megaleporinus muyscorum</em></a> by Grupo de Ictiología de la Universidad de Antioquia, licensed under <a href='https://creativecommons.org/licenses/by/4.0/' target='_blank'>CC by 4.0</a>. The material was not modified</li><li>Atherinopsidae: <a href='https://www.phylopic.org/images/4e31c12d-1c23-4a64-bcc6-0730ab6d0617/odontesthes' target='_blank'><em>Odontesthes</em></a> by Lily Hughes, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Bagridae: <a href='https://www.phylopic.org/images/9d7ceefe-f1e3-4f1d-bec2-0a52271ac236/bagridae' target='_blank'><em>Bagridae</em></a> by Cee Nell, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Centrarchidae: <a href='https://www.phylopic.org/images/ccee99b3-052a-4e2d-9db2-80a3a75484f0/pomoxis-nigromaculatus' target='_blank'><em>Pomoxis nigromaculatus</em></a> by Rene Martin, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Centropomidae: <a href='https://www.phylopic.org/images/e0cd6df8-23c0-4973-aa3a-8ec042a4cea2/centropomidae' target='_blank'><em>Centropomidae</em></a> by Cee Nell, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Channidae: <a href='https://www.phylopic.org/images/8e3b6948-6337-4802-9d56-067ce1f41a53/channa-argus' target='_blank'><em>Channa argus</em></a> by Christopher Kenaley, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Characidae: <a href='https://www.phylopic.org/images/a3cd77e8-f869-4506-90de-38c45daefb6d/hemibrycon-boquiae' target='_blank'><em>Hemibrycon boquiae</em></a> by Grupo de Ictiología de la Universidad de Antioquia, licensed under <a href='https://creativecommons.org/licenses/by/4.0/' target='_blank'>CC by 4.0</a>. The material was not modified</li><li>Cichlidae: <a href='https://www.phylopic.org/images/e7e9a112-9daa-4d72-be90-bc0bdd5473f6/cichlidae' target='_blank'><em>Cichlidae</em></a> by Cee Nell, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Cyprinidae: <a href='https://www.phylopic.org/images/25a76141-cca2-4b56-bd9f-9f5fdcd8efa1/cyprinus-carpio' target='_blank'><em>Cyprinus carpio</em></a> by Carlos Cano-Barbacil, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Erythrinidae: <a href='https://www.phylopic.org/images/42887a53-f47d-4d00-8251-23b902bd8608/hoplias-malabaricus' target='_blank'><em>Hoplias malabaricus</em></a> by Grupo de Ictiología de la Universidad de Antioquia, licensed under <a href='https://creativecommons.org/licenses/by/4.0/' target='_blank'>CC by 4.0</a>. The material was not modified</li><li>Esocidae: <a href='https://www.phylopic.org/images/e8173d95-0261-4984-a5c2-5111f213e1d4/esox-lucius' target='_blank'><em>Esox lucius</em></a> by Russell Engelman, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Ictaluridae: <a href='https://www.phylopic.org/images/4317e447-87bb-42b9-94bc-8c910afdbe4b/ictalurus-punctatus' target='_blank'><em>Ictalurus punctatus</em></a> by Carlos Cano-Barbacil, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Latidae: <a href='https://www.phylopic.org/images/6821092f-e3f0-46ad-970c-a1969cc9e5d6/lates-niloticus' target='_blank'><em>Lates niloticus</em></a> by Jessica Rick, licensed under <a href='https://creativecommons.org/licenses/by-nc/3.0/' target='_blank'>CC BY-NC 3.0</a>. The material was not modified</li><li>Lotidae: <a href='https://www.phylopic.org/images/9e87f9e3-a1a6-49e0-9037-924d81ab53a3/lota-lota' target='_blank'><em>Lota lota</em></a> by Ellen Edmonson (illustration) and Timothy J. Bartley (silhouette), licensed under <a href='https://creativecommons.org/licenses/by-sa/3.0/' target='_blank'>CC BY-SA 3.0</a>. The material was not modifed</li><li>Osmeridae: <a href='https://www.phylopic.org/images/cd4174a2-9c45-47aa-b4db-b88f9eece5e4/osmerus-eperlanus' target='_blank'><em>Osmerus eperlanus</em></a> by Fidji Berio, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Osteichthyes: <a href='https://www.phylopic.org/images/66674c55-bafe-4dae-ae9a-3f84d7d7fd8a/squalius-cephalus' target='_blank'><em>Squalius cephalus</em></a> by S. Ginot, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Percichthyidae: <a href='https://www.phylopic.org/images/3210ebc7-32cb-463e-9370-3be415be9e02/percichthyidae' target='_blank'><em>Percichthyidae</em></a> by Cee Nell, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Percidae: <a href='https://www.phylopic.org/images/b1649388-f901-41b2-9c31-0ed9b1057cd9/sander-vitreus' target='_blank'><em>Sander vitreus</em></a> by NOAA Great Lakes Environmental Research Laboratory (illustration) and Timothy J. Bartley (silhouette), licensed under <a href='https://creativecommons.org/licenses/by-sa/3.0/' target='_blank'>CC BY-SA 3.0</a>. The material was not modified</li><li>Pimelodidae: <a href='https://www.phylopic.org/images/2999e387-471d-490c-b838-fb15a19dfe91/phractocephalus-hemioliopterus' target='_blank'><em>Phractocephalus hemioliopterus</em></a> by Andy Wilson, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Plecoglossidae: <a  href='https://www.phylopic.org/images/0cdeb298-61d5-472f-b255-953480a4d8c8/plecoglossus-altivelis' target='_blank'><em>Plecoglossus altivelis</em></a> by Cee Nell, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Salmonidae: <a href='https://www.phylopic.org/images/5ffc1b34-6602-4f2a-a9e1-9078d7c43627/oncorhynchus-gorbuscha' target='_blank'><em>Oncorhynchus gorbuscha</em></a> by Tessa Rehill, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Serrasalmidae: <a  href='https://www.phylopic.org/images/173896a9-b0ee-4179-9056-44165c07ed63/metynnis' target='_blank'><em>Metynnis</em></a> by Andy Wilson, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Siluridae: <a href='https://www.phylopic.org/images/7cb387d7-9b23-4e5c-a007-2d938339cb8a/silurus-glanis' target='_blank'><em>Silurus glanis</em></a> by Carlos Cano-Barbacil, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li><li>Sinipercidae: <a  href='https://www.phylopic.org/images/95a6de68-d304-48ea-9235-c8b3661ae4c1/sinipercidae' target='_blank'><em>Sinipercidae</em></a> by Cee Nell, licensed under <a href='https://creativecommons.org/publicdomain/zero/1.0/' target='_blank'>CC0 1.0</a></li></ul></span>"
                        },
                    ],
                    activeOnLoad: false
                },
            ]
        },        
        FishAsFoodSankey: {
            paragraph1: "Explore inland recreational <span class='tooltip-group'><span class='tooltip-span'>consumable harvest</span><span id='consumable-tooltip' class='tooltiptext'>Calculated from species-specific estimates of literature-based edible portion (%) of a given fish</span></span> for the five families with the largest recreational consumable harvest: <span class='scientificName'>Cyprinidae</span> (minnows and carps), <span class='scientificName'>Percidae</span> (perch), <span class='scientificName'>Salmonidae</span> (salmon, trout, grayling, and whitefish), <span class='scientificName'>Bagridae</span> (bagrid catfish), and <span class='scientificName'>Centrarchidae</span> (sunfishes). Total harvest is split by family (at left), species/co-reported taxa (in the center), and country (at right). Hover over the chart to see the harvest totals in kilograms. Use the toggles to include or exclude China from the chart and to sort the chart from highest to lowest harvest totals.",
            toggleData: {
                showChina: {
                    value: true,
                    color: 'var(--usgs-blue)',
                    label: "include China"
                },
                sortHarvest: {
                    value: false,
                    color: 'var(--usgs-blue)',
                    label: "sort by harvest totals"
                }
            }
        },
        FishAsFoodClimate: {
            paragraph1: `Inland recreational fishing is often considered a leisure-driven activity, yet <a href='/${viteAppAssetURl}/earth-in-flux/rec-fish-as-food/rec-fish-consumption' target='_blank'>recreationally harvested fish</a> can be an important food source. Depending on what species are fished, this consumption may be more or less vulnerable to climate change. Use the chart below to explore the global implications of inland recreational fishing as a source of food in a changing climate.`,
            image: "fish_as_food_continent_map.png",
            imageAlt: "Map of the continents for which data are shown",
            defaultInfo: {
                image: "fish_as_food_continent_map.png",
                text: "In the chart below, each circle represents a country. The circles are colored by continent, and the size of each circle is proportional to the number of recreational fishers in each country.",
                prompt_desktop: "Click on a circle to learn more about the climate vulnerability of recreational fish consumption in each country. Click outside the circle to return to the original view.",
                prompt_mobile: "Tap on a circle to learn more about the climate vulnerability of recreational fish consumption in each country.",
                name: null,
                continent: null
            }
        },
        WildfireAerosols: {
            paragraph1: `Each layer of the <a href='/${viteAppAssetURl}/earth-in-flux/fire-in-ice/glacier-scan' target='_blank'>collected snow core</a> contains more than just packed snow. Particulates from the air, like dust, deposit on the surface of the glacier, along with tiny airborne particles called aerosols. Over time, the deposited particulates and aerosols are preserved in the glacier. If the aerosols are dark in color, the glacier absorbs more heat and melts more quickly. These dark aerosols include black carbon, or soot, that is generated when vehicles and industrial activities burn fossil fuels or when wildfires burn vegetation.`,   
            paragraph2: "Can we tell if any of the dark aerosols in the snow core came from wildfires? While black carbon does not have a chemical signature, three sugars—mannnosan, galactosan, and levoglucosan—are only produced when vegetation burns. While these sugars are not dark aerosols themselves, they travel in smoke plumes with the dark aerosols and are deposited alongside them. These sugars are present throughout the core, which tells us that some of the deposited material in the snow came from wildfires.",   
            paragraph3: `Scientists use the ratio of levoglucosan to the sum of mannosan and galactosan to distinguish between types of vegetation that burned. Alaska's forests are dominated by <span class='tooltip-group'><span class='tooltip-span'>softwoods</span><span id='softwoods-tooltip' class='tooltiptext'>Conifers, like pine, spruce, and firs.</span></span>, and <a href='/${viteAppAssetURl}/earth-in-flux/fire-in-ice/regional-fires' target='_blank'>regional fires</a> likely deposit aerosols on the Juneau Icefield that are captured in the core.`,
            paragraph4: "However, there are also markers of <span class='tooltip-group'><span class='tooltip-span'>hardwood</span><span id='hardwoods-tooltip' class='tooltiptext'>Broadleaved trees, like ash, oak, birch, beech, alder, and teak.</span></span> combustion, which suggests that wildfire aerosols are transported to the icefield from much farther afield. One possible source is wildfires from East Asia, where hardwood forests are more abundant.",
            heading: "What am I looking at?",
            explanation1: `As scientists <a href='/${viteAppAssetURl}/earth-in-flux/fire-in-ice/glacier-scan' target='_blank'>collected the snow core</a>, they carefully stored each ten-centimeter section for transport off the Juneau Icefield. The retrieved samples were analyzed in a laboratory for particle counts, major ions, stable isotopes of oxygen and hydrogen, and the three sugars that are markers of biomass combustion—mannosan, galactosan, and levoglucosan.`,
            explanation2: "In this visual representation of the core, the darker grey shows layers of the snow that had high amounts of particulate matter. The high-particulate layer present in both 2015 and 2016 likely represents the summer melt surface."
        },
        FindexThreatDefinitions: {
            paragraph1: "Inland fisheries face a range of threats driven by human activities and environmental change. These pressures can affect aquatic habitats, water quality, and fish populations. Explore five categories of threats to learn what drives them and how they affect inland fisheries.",
            iconData: [
                {
                    iconTitle: "Habitat",
                    iconID: "habitat",
                    icon: "habitat-29572",
                    iconText: "Fish need healthy, connected waters to spawn, feed, and grow. Dams, drained wetlands, cleared forests, and heavy water use reshape rivers and lakes, cutting fish off from the places they depend on.",
                    subThreatData: [
                        {
                            subThreat: "Dams",
                            subThreatIcon: "H-dams-3916913",
                            subThreatText: "Hydropower, water supply for growing urban populations, flood control for urban and agricultural areas, and irrigation for agriculture."
                        },
                        {
                            subThreat: "Wetland drainage",
                            subThreatIcon: "H-wetland-drainage-34007",
                            subThreatText: "Wetlands are drained directly for agriculture, aquaculture, housing, industry, and recreation, and indirectly when nearby water withdrawals lower the water table. Climate-driven evaporation compounds the loss. Upstream drainage also reduces flows to downstream wetlands, raising the risk of salinization near coasts."
                        },
                        {
                            subThreat: "Deforestation and associated runoff",
                            subThreatIcon: "H-deforestation-4811548",
                            subThreatText: "Demand for timber products and conversion of forests to agricultural and urban land."
                        },
                        {
                            subThreat: "Riparian degradation",
                            subThreatIcon: "H-riparian-6291",
                            subThreatText: "Land use changes such as deforestation, watershed regulation, hydropower, agriculture, irrigation, water withdrawals, urbanization, and pollution."
                        },
                        {
                            subThreat: "Agricultural water use",
                            subThreatIcon: "H-agriculture-5491228",
                            subThreatText: "Growing demand for food, fiber, feed, and fuel. Local water supply can be reduced through irrigation, affecting fish habitats."
                        },
                        {
                            subThreat: "Urban water use",
                            subThreatIcon: "H-urban-7019973",
                            subThreatText: "Population growth, economic development, wastewater intake and discharge, recreational fishing, and river transport alter natural flows. Wastewater returning to local waterbodies can add to the impact."
                        },
                        {
                            subThreat: "Industrial water use",
                            subThreatIcon: "H-industrial-659077",
                            subThreatText: "Economic development and industrial wastewater can strain local water resources."
                        }
                    ]
                },
                {
                    iconTitle: "Pollution",
                    iconID: "pollution",
                    icon: "pollution-6278",
                    iconText: "What happens on land doesn't stay on land. Runoff and wastewater from farms, cities, industry, aquaculture, and mining carry nutrients and contaminants into rivers and lakes, where they can harm fish and the people who eat them.",
                    subThreatData: [
                        {
                            subThreat: "Agricultural effluents",
                            subThreatIcon: "P-agriculture-3634692",
                            subThreatText: "Agricultural runoff can contain elevated nutrients, affecting downstream water resources."
                        },
                        {
                            subThreat: "Urban wastewater",
                            subThreatIcon: "P-urban-waste-25100",
                            subThreatText: "Urbanizing populations around the world are increasing municipal wastewater discharge."
                        },
                        {
                            subThreat: "Industrial effluents",
                            subThreatIcon: "P-industrial-200223",
                            subThreatText: "Industrial runoff and point-source pollution can affect local and downstream waterbodies."
                        },
                        {
                            subThreat: "Aquaculture effluents",
                            subThreatIcon: "P-aquaculture-23815",
                            subThreatText: "Expanding aquaculture production as fish supply shifts from capture fisheries to farming."
                        },
                        {
                            subThreat: "Pharmaceuticals",
                            subThreatIcon: "P-pharmaceuticals-164439",
                            subThreatText: "Growing use of pharmaceuticals in human and veterinary medicine."
                        },
                        {
                            subThreat: "Oil or gas exploration",
                            subThreatIcon: "P-oil-gas-3723255",
                            subThreatText: "Water used for oil or gas extraction can contain various water quality contaminants."
                        },
                        {
                            subThreat: "Plastics",
                            subThreatIcon: "P-plastics-2484814",
                            subThreatText: "Increasing production and consumption of plastics."
                        },
                        {
                            subThreat: "Mining",
                            subThreatIcon: "P-mining-7248907",
                            subThreatText: "Water used for mining exploration and extraction can affect downstream systems, and direct input of oil or gas from leaks can directly impact water quality."
                        }
                    ]
                },
                {
                    iconTitle: "Climate and weather",
                    iconID: "climate-and-weather",
                    icon: "climate-6795675",
                    iconText: "As the climate changes, so do the waters fish call home. Warmer temperatures, more extreme droughts and floods, shifting winds, and shorter ice seasons are changing where and how fish can thrive.",
                    subThreatData: [
                        {
                            subThreat: "Change in water temperature",
                            subThreatIcon: "CW-temperature-1979336",
                            subThreatText: "Shifting weather patterns, snowmelt quantities, and warming temperatures affect local waterbody temperature."
                        },
                        {
                            subThreat: "Drought",
                            subThreatIcon: "CW-drought-002D5E",
                            subThreatText: "More frequent and intense periods of abnormally low precipitation or reduced snowmelt."
                        },
                        {
                            subThreat: "Change in flooding",
                            subThreatIcon: "CW-flooding-1760091",
                            subThreatText: "Changes in the frequency and intensity of precipitation events and sea-level rise."
                        },
                        {
                            subThreat: "Change in wind patterns",
                            subThreatIcon: "CW-wind-7479181",
                            subThreatText: "Warmer water may reduce land–water temperature contrasts, increasing the variability or severity of winds and changing lake stratification processes and productivity."
                        },
                        {
                            subThreat: "Change in ice cover",
                            subThreatIcon: "CW-ice-cover-5314745",
                            subThreatText: "Rising temperatures, changes in season timing, and increased rain-on-snow events affect ice cover on lakes and streams."
                        }
                    ]
                },
                {
                    iconTitle: "Invasive species",
                    iconID: "invasive-species",
                    icon: "invasive-species-973991",
                    iconText: "When non-native species take hold in new waters, native fish can lose out. Invaders compete for food, prey on native species, spread disease, and disrupt entire food webs.",
                    subThreatData: [
                        {
                            subThreat: "Invasive non-native species",
                            subThreatIcon: "invasive-species-973991",
                            subThreatText: "Species spread accidentally on boats, equipment, and stocked fish, or are introduced deliberately for aquaculture, stocking, or biological control."
                        }
                    ]
                },
                {
                    iconTitle: "Fishing pressure",
                    iconID: "fishing-pressure",
                    icon: "fishing-pressure-4847095",
                    iconText: "Fishing feeds and supports millions of people, but when fish are caught faster than populations can recover, stocks decline. Commercial, subsistence, and recreational fishing all play a role.",
                    subThreatData: [
                        {
                            subThreat: "Overfishing",
                            subThreatIcon: "fishing-pressure-4847095",
                            subThreatText: "Rising demand from rural populations, food-insecure communities, and commercial markets. Many inland fisheries are open access, increasing the risk of heavy fishing and illegal gear use. Recreational fishing can also heavily target large species."
                        }
                    ]
                }
            ]
        },
        FindexThreatSankey: {
            paragraph1: 'Land use change is the biggest threat to inland fisheries, with dams, wetland drainage, and deforestation and associated changes in runoff as the largest contributors. View the Sankey plot below to see the biggest threats to inland fisheries across threat categories of habitat, pollution, climate and weather, invasive species, and fishing pressure.'
        },
        FindexGlobalThreats: {
            heading1:"What are the top threats to inland fisheries?",
            paragraph1: `While threats related to habitat are the <a href='/${viteAppAssetURl}/earth-in-flux/findex/inland-fish-total-threats' target='_blank'>overall highest threat</a> to inland fisheries worldwide, the top threat for each river basin varies globally. The map below shows the top threat to inland fisheries in each river basin.`,
            mapData: {
                habitat: {
                    visible: true,
                    path: 'Habitat_threat_by_basin.png',
                    color: 'var(--color-habitat-top)',
                    order: 1,
                    label: "Habitat"
                },
                pollution: {
                    visible: true,
                    path: 'Pollution_threat_by_basin.png',
                    color: 'var(--color-pollution-top)',
                    order: 2,
                    label: "Pollution"
                },
                climate_and_weather: {
                    visible: true,
                    path: 'Climate_and_weather_threat_by_basin.png',
                    color: 'var(--color-climate-and-weather-top)',
                    order: 3,
                    label: "Climate and weather"
                },
                invasive_species: {
                    visible: true,
                    path: 'Invasive_species_threat_by_basin.png',
                    color: 'var(--color-invasive-species-top)',
                    order: 4,
                    label: "Invasive species"
                },
                fishing_pressure: {
                    visible: true,
                    path: 'Fishing_pressure_threat_by_basin.png',
                    color: 'var(--color-fishing-pressure-top)',
                    order: 5,
                    label: "Fishing pressure"
                },
            },
            heading2:"Where are threats higher and lower?",
            paragraph2: "Similarly, the severity of threats varies globally. Explore the tabs below to learn more.",
            subThreatHeading1: "Drivers",
            subThreatHeading2: "Potential impacts",
            subThreatHeading3: "Sectors involved",
            subThreatHeading4: "Mitigation and adaptation opportunities",
            subThreatHeading5: "How the threat is measured",
            tabData: [
                {
                    tabTitle: "Habitat",
                    tabContentTitle: "Habitat",
                    tabContentTitleID: "habitat",
                    tabIcon: "habitat-29572",
                    tabText: "Habitat loss and alteration is the greatest threat to inland fisheries worldwide. Dams, wetland drainage, deforestation, riparian degradation, and water withdrawals for agriculture, cities, and industry change how water, sediment, and nutrients move through river basins, cutting fish off from the spawning, nursery, and feeding areas they depend on. Select a subthreat category above to learn more.",
                    tabMapImageAlt: "Habitat alt.",
                    tabLegendImageAlt: "Legend",
                    tabColor: "var(--color-habitat)",
                    subThreatPrefix: "H",
                    subThreatData: [
                        {
                            subThreat: "Dams",
                            subThreatIcon: "H-dams-3916913",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Hydropower, water supply for growing urban populations, flood control for urban and agricultural areas, and irrigation for agriculture.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Dams block river connectivity, interrupting sediment and nutrient transport and disrupting fish migration and recruitment. Productivity declines upstream and downstream, floodplains flood less often, and low flows increase. Reservoirs replace river fish with lake-adapted species, and drawdown and refilling destroy aquatic vegetation and breeding areas. Altered flows also remove thermal cues, reducing fish growth and reproduction. These effects are well documented (<a href='https://archive.org/details/impoundedriversp0000pett/mode/2up' target='_blank'>Petts 1984</a>), including in large tropical rivers (<a href='https://doi.org/10.1126/science.aac7082' target='_blank'>Winemiller et al. 2016</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Energy, water supply, and flood mitigation sectors.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Build fish passage facilities; work with reservoir managers to ensure minimum flows and timely flood releases; stock species suited to reservoirs, or stock young fish where recruitment is limited; and prepare fishing communities for changing fisheries. Options for reducing dam impacts are reviewed by <a href='https://doi.org/10.1016/j.watbs.2022.100009' target='_blank'>Twardek et al. (2022)</a>.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Average of three dam impact metrics per stream kilometer: degree of fragmentation (how much dams divide the river network), degree of regulation (reservoir storage relative to annual river flow), and sediment trapping (share of the natural sediment load held behind dams) (<a href='https://doi.org/10.6084/m9.figshare.7688801' target='_blank'>Grill and Lehner 2019</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Wetland drainage",
                            subThreatIcon: "H-wetland-drainage-34007",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Wetlands are drained directly for agriculture, aquaculture, housing, industry, and recreation, and indirectly when nearby water withdrawals lower the water table. Weak legislation and governance enable drainage, and climate-driven evaporation compounds it. Upstream drainage also reduces flows to downstream wetlands, raising the risk of salinization near coasts.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Loss of habitat, especially floodplain wetlands used for spawning, reduces species diversity; endemic species in small wetlands are most at risk. Shallower water lowers water quality and increases temperature extremes, and impacts can extend downstream. Consolidating small wetlands into one larger, deeper wetland may support more fish but provides poorer habitat for invertebrates, birds, and submerged plants (<a href='https://doi.org/10.1890/ES14-00494.1' target='_blank'>McCauley et al. 2015</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Agricultural, industrial, municipal, and water sectors.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Restore wetlands, recognizing that restored wetlands may not regain full ecological function (<a href='https://doi.org/10.1371/journal.pbio.1001247' target='_blank'>Moreno-Mateos et al. 2012</a>). Preventing drainage of carbon-rich wetlands such as peatlands also reduces greenhouse gas emissions (<a href='https://doi.org/10.1038/s41558-020-0738-8' target='_blank'>Goldstein et al. 2020</a>). Restoration can pay off: re-flooding the Waza Logone floodplain in Cameroon was estimated to add $2.5 million per year to the regional economy, including through restored fisheries (<a href='https://doi.org/10.2305/IUCN.CH.2004.WANI.3.en' target='_blank'>Emerton and Bos 2004</a>).",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Percent of basin area where wetlands were lost between 1700 and 2020 (<a href='https://doi.org/10.5281/zenodo.7616651' target='_blank'>Fluet-Chouinard et al. 2023</a>). Basins with no net wetland loss score zero.",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Deforestation and associated runoff",
                            subThreatIcon: "H-deforestation-4811548",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Demand for timber products and conversion of forests to agricultural and urban land.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Cleared land erodes easily, sending sediment into waterways and degrading habitat. Sedimentation reduces food production and invertebrate densities, abrades gills, lowers spawning success for species that spawn on gravel and rock, and impairs foraging by reducing visibility (<a href='https://doi.org/10.1002/hyp.7940' target='_blank'>Kemp et al. 2011</a>). Loss of shade warms water, stressing or killing temperature-sensitive fish (<a href='https://doi.org/10.1371/journal.pone.0196560' target='_blank'>Ilha et al. 2018</a>), and reduces the woody debris some fish rely on for cover and food (<a href='https://doi.org/10.1002/aqc.3270050206' target='_blank'>Gurnell et al. 1995</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Forestry, agriculture, municipal planners, and governments.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Limit deforestation in sensitive or ecologically valuable basins, protect riparian buffers, and use logging techniques that minimize erosion. Forested riparian buffers have been shown to reduce the effects of deforestation on fish communities (<a href='https://doi.org/10.1890/08-0050.1' target='_blank'>Lorion and Kennedy 2009</a>).",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Percent of basin area where forest cover was lost since 2000 (<a href='https://doi.org/10.1126/science.1244693' target='_blank'>Hansen et al. 2013</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Riparian degradation",
                            subThreatIcon: "H-riparian-6291",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Land use changes such as deforestation, watershed regulation, hydropower, agriculture, irrigation, water withdrawals, urbanization, and pollution.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Fish responses vary with the type of water body and the nature and degree of disturbance (<a href='https://doi.org/10.1038/s41467-020-14624-2' target='_blank'>Kao et al. 2020</a>). Opportunistic species may thrive while specialized species decline, but almost any change on land alters the composition and structure of fish communities. In Lake Victoria, land degradation and species introductions have sharply reduced fish diversity even as overall productivity increased (<a href='https://doi.org/10.1079/9781845934149.0309' target='_blank'>Kolding et al. 2008</a>; <a href='https://doi.org/10.1139/cjfas-2015-0130' target='_blank'>van Zwieten et al. 2016</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Forestry, agriculture, energy, industry, municipal planners, and governments.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Recognition among policymakers and planners that changes to riparian land inevitably affect waterways. Planned degradation can be partly offset through technical interventions, and decades of restoration experience show what works: progress depends on willingness, environmental legislation, and economic trade-offs.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Percent of land within riparian buffers that is built-up or cropland, based on 2019 land cover data (<a href='https://doi.org/10.5281/zenodo.3939050' target='_blank'>Buchhorn et al. 2020</a>). Buffer widths scale with stream and lake size.",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Agricultural water use",
                            subThreatIcon: "H-agriculture-5491228",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Growing demand for food, fiber, feed, and fuel. Agricultural markets are global, but water is managed locally, and irrigation water is rarely managed to maintain downstream connectivity. Large transboundary water resources are often governed by outdated, inflexible rules.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Withdrawals reduce or destabilize downstream flows, and the associated infrastructure, such as dams and reservoirs, can cut off floodplains and block fish movement and migration (<a href='https://www.fao.org/3/CA2675EN/ca2675en.pdf' target='_blank'>Gregory et al. 2018</a>). This can reduce biodiversity, alter fish communities, increase fish mortality, and reduce catches, contributing to malnutrition (<a href='https://doi.org/10.1080/00207230600963296' target='_blank'>Degeorges and Reilly 2006</a>). However, irrigation reservoirs can also moderate flow extremes and create fish habitat where fisheries would not otherwise exist (<a href='https://www.fao.org/3/CA2675EN/ca2675en.pdf' target='_blank'>Gregory et al. 2018</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Farmers, ranchers, water managers, government officials, local and downstream communities, and environmental professionals.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Redesign and operate irrigation delivery and storage infrastructure to improve connectivity and flows, create or improve fish habitat and refuges in and around irrigated systems, and revise policies and regulations to support these changes (<a href='https://www.fao.org/3/CA2675EN/ca2675en.pdf' target='_blank'>Gregory et al. 2018</a>). Raising fish in irrigation reservoirs and rice fields can also help offset fisheries lost to altered flows (<a href='https://www.fao.org/4/a0444e/A0444E09.htm#ch7' target='_blank'>Peterson and Kalende 2006</a>; <a href='https://digitalarchive.worldfishcenter.org/server/api/core/bitstreams/e2148d94-7037-4b51-b162-73f9705f137e/content' target='_blank'>Baran 2005</a>).",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Water consumed for irrigation and livestock in 2020, averaged across five climate models, relative to river discharge in the basin (<a href='https://doi.org/10.7910/dvn/viqeab' target='_blank'>Khan et al. 2022</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Urban water use",
                            subThreatIcon: "H-urban-7019973",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Population growth, economic development, wastewater intake and discharge, recreational fishing, and river transport alter natural flows. Weak legislation on treating water before it returns to rivers adds to the problem.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Discharges from water treatment add pollution and excess nutrients to rivers. Withdrawals and channel modifications for transport alter flows and cut fish off from floodplains, changing species composition and reducing access to spawning areas. Diverting river water to fill artificial recreational fishing lakes further disconnects rivers from their floodplains.",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Water management and water-use sectors, governments, and commercial water users.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Enforce wastewater treatment before release into rivers, maintain environmental flows in urban areas (particularly around transport infrastructure), and restore urban rivers and fish refuges. The <a href='https://www.wildtrout.org/assets/img/general/Urban-River-Toolkit.pdf' target='_blank'>Trout in the Town Urban River Toolkit</a> offers practical restoration measures.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Water consumed for domestic use and electricity generation in 2020, averaged across five climate models, relative to river discharge in the basin (<a href='https://doi.org/10.7910/dvn/viqeab' target='_blank'>Khan et al. 2022</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Industrial water use",
                            subThreatIcon: "H-industrial-659077",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Economic development combined with weak environmental policy and legislation, and poor coordination between water and fisheries management authorities.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Withdrawals reduce environmental flows, increasing habitat variability and disrupting spawning routes. Fish species composition shifts, and fisheries see smaller, less predictable catches of different species.",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Water managers, governments, and upstream and downstream communities and authorities.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Give greater weight to environmental flows through integrated water resources management and establish fish refuges near heavy industrial areas.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Water consumed for mining and manufacturing in 2020, averaged across five climate models, relative to river discharge in the basin (<a href='https://doi.org/10.7910/dvn/viqeab' target='_blank'>Khan et al. 2022</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        }
                    ]
                },
                {
                    tabTitle: "Pollution",
                    tabContentTitle: "Pollution",
                    tabContentTitleID: "pollution",
                    tabIcon: "pollution-6278",
                    tabText: "Pollution enters inland waters from farms, cities, industry, aquaculture, and resource extraction. Excess nutrients deplete oxygen, while pesticides, heavy metals, pharmaceuticals, petroleum, and plastics can poison fish, impair their growth, behavior, and reproduction, and contaminate the fish people eat. Select a subthreat category above to learn more.",
                    tabMapImageAlt: "Pollution alt.",
                    tabLegendImageAlt: "Legend",
                    tabColor: "var(--color-pollution)",
                    subThreatPrefix: "P",
                    subThreatData: [
                        {
                            subThreat: "Agricultural effluents",
                            subThreatIcon: "P-agriculture-3634692",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Food demand, economic policies, and weak environmental regulation and enforcement.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Pesticides, including organochlorines, organophosphates, and pyrethroid insecticides, are highly toxic to fish. They cause oxidative stress, alter blood chemistry, suppress immune function, damage organs such as the liver and kidneys, and change fish behavior. These effects have been documented in Nile tilapia (<a href='https://doi.org/10.1007/s11259-022-09961-0' target='_blank'>Abdelmagid et al. 2022</a>), and the toxicity of the herbicide glyphosate to aquatic organisms is reviewed by <a href='https://ph03.tci-thaijo.org/index.php/ahstr/article/view/1523' target='_blank'>Thanomsit et al. (2020)</a>.",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Agriculture, governments, rural communities, and academia.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Catchment management; better regulation, monitoring, and enforcement; and improved biological information.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Agricultural grey water footprint: the volume of fresh water needed to dilute nitrogen from crop fertilizers to meet water quality standards (<a href='https://doi.org/10.1073/pnas.1109936109' target='_blank'>Hoekstra and Mekonnen 2012</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Urban wastewater",
                            subThreatIcon: "P-urban-waste-25100",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Urbanizing populations around the world are increasing municipal wastewater discharge.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Urban wastewater can contain pollutants harmful to fish. Organic material and nutrients cause eutrophication and lower dissolved oxygen; toxic contaminants such as ammonia, chloride, and industrial chemicals can be acutely toxic; and hormones, pharmaceuticals, and personal care products can disrupt endocrine systems and impair reproduction (<a href='https://doi.org/10.1002/ieam.1364' target='_blank'>Tetreault et al. 2013</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Municipal governments and regulatory agencies.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Expand effective wastewater treatment and strengthen regulation of municipal wastewater. Improved wastewater management can substantially reduce the threat to downstream fish communities (<a href='https://doi.org/10.1139/f2011-096' target='_blank'>Holeton et al. 2011</a>).",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Domestic grey water footprint: the volume of fresh water needed to dilute household wastewater pollutants to meet water quality standards (<a href='https://doi.org/10.1073/pnas.1109936109' target='_blank'>Hoekstra and Mekonnen 2012</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Industrial effluents",
                            subThreatIcon: "P-industrial-200223",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Economic development combined with weak environmental policy and legislation.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Industrial effluents may raise nutrient, toxin, and salt concentrations, depleting oxygen and altering food webs. Fish can suffer acute toxicity and reduced spawning success, and fisheries lose value as products become contaminated.",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Governments (environment and industry ministries), fisheries, consumers, and downstream communities.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Stricter, industry-specific standards for treating effluents and capacity building on their environmental and socioeconomic effects. Communities can adapt through fish consumption guidelines and innovative or nature-based solutions for water purification.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Industrial grey water footprint: the volume of fresh water needed to dilute industrial pollutants to meet water quality standards (<a href='https://doi.org/10.1073/pnas.1109936109' target='_blank'>Hoekstra and Mekonnen 2012</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Aquaculture effluents",
                            subThreatIcon: "P-aquaculture-23815",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Expanding aquaculture production as fish supply shifts from capture fisheries to farming.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Effluents can contain nutrients, solids, pharmaceuticals, pathogens, and escaped farmed organisms. Nutrients drive eutrophication, and solids increase turbidity and sedimentation. Most aquaculture pharmaceuticals are nontoxic for fish, but some, such as hormones used for sex reversal, have sublethal effects. Pathogens can severely affect wild fish, especially when non-native pathogens are introduced or farmed fish act as reservoirs for native ones.",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Aquaculture, commercial and recreational fishing industries, and government resource managers.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Effluents from tank systems and, to a lesser extent, pond systems can be treated. Cage and pen systems are open to the environment and offer no such opportunity.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Combined nitrogen and phosphorus loads from aquaculture, estimated with the IMAGE-GNM global nutrient model (<a href='https://doi.org/10.17026/dans-zgs-9k9m' target='_blank'>Beusen and PBL 2016</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Pharmaceuticals",
                            subThreatIcon: "P-pharmaceuticals-164439",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Growing use of pharmaceuticals in human and veterinary medicine.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Depending on the species and the type and concentration of pharmaceutical, fish behavior, reproduction, development, and genetics may be harmed (<a href='https://doi.org/10.1002/ieam.1364' target='_blank'>Tetreault et al. 2013</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Pharmaceutical industry, medical sector, livestock farmers, and health and environmental agencies.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Reduce pharmaceutical use in animals, ensure proper disposal, and improve detection in waterways.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Modeled river concentrations of diclofenac, a widely used anti-inflammatory drug, as an indicator of pharmaceutical contamination (<a href='https://doi.org/10.1016/j.envint.2020.105993' target='_blank'>Acuña et al. 2020</a>). Lakes are not included.",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Oil or gas exploration",
                            subThreatIcon: "P-oil-gas-3723255",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Weak legislation, limited enforcement, and lack of stakeholder involvement.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Polycyclic aromatic hydrocarbons (PAHs) from petroleum disrupt salmonid embryo development, impair swimming, and cause malformations. Fish exposed to oil or gas extraction effluents or direct leaks become more susceptible to parasites and to eye and heart problems. Water bodies affected by crude oil show marked declines in fish diversity (<a href='https://doi.org/10.1590/1982-0224-2021-0001' target='_blank'>Azevedo-Santos et al. 2021</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Oil industry, international corporations, and governments.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Better regulation, monitoring, and enforcement, along with stakeholder consultation.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Percent of basin area overlying oil and gas fields (<a href='https://doi.org/10.18141/1427300' target='_blank'>Sabbatino 2018</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Plastics",
                            subThreatIcon: "P-plastics-2484814",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Increasing production and consumption of plastics.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Plastics are found in waters throughout the world. As they break down into microplastics (pieces smaller than 5 millimeters), they spread through the water column and sediments and are ingested by invertebrates, fish, birds, and mammals (<a href='https://doi.org/10.1016/j.envpol.2020.115083' target='_blank'>Spanjer et al. 2020</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Manufacturers, industry, and environmental agencies.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Reduce plastic use and remove plastics and microplastics from waterways.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Mismanaged plastic waste entering rivers, estimated at 40,760 river mouths worldwide and distributed upstream in proportion to human population (<a href='https://doi.org/10.6084/m9.figshare.4725541.v6' target='_blank'>Lebreton and Reisser 2020</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Mining",
                            subThreatIcon: "P-mining-7248907",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Poverty and lack of formal employment, weak legislation, and limited enforcement.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Mining increases erosion and siltation, and toxic sludge and chemicals used to separate minerals can reach rivers. Heavy metals can impair fish growth, reproduction, and nervous system function, and cyanide can cause delayed mortality, impaired swimming, disrupted respiration, and altered growth (<a href='https://pubs.usgs.gov/publication/5200123' target='_blank'>Eisler 1991</a>). In southeastern Brazil, the collapse of mine tailings dams sent toxic waste into the Doce River, severely affecting fish diversity (<a href='https://doi.org/10.1016/j.pecon.2017.06.002' target='_blank'>Carmo et al. 2017</a>; <a href='https://doi.org/10.1016/j.ecoenv.2020.110165' target='_blank'>Weber et al. 2020</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Industrial and artisanal mining, governments, and environmental agencies.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Better regulation, monitoring, and control of tailings disposal facilities.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Percent of basin area covered by mining operations, including large-scale, artisanal, and small-scale mines (<a href='https://doi.org/10.1594/pangaea.942325' target='_blank'>Maus et al. 2022</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        }
                    ]
                },
                {
                    tabTitle: "Climate and weather",
                    tabContentTitle: "Climate and weather",
                    tabContentTitleID: "climate-and-weather",
                    tabIcon: "climate-6795675",
                    tabText: "Climate change is altering the conditions inland fish depend on. Warmer water, more intense droughts and floods, shifting wind patterns, and shorter periods of lake ice cover affect fish physiology, reproduction, and habitat, and are changing where species can survive. Select a subthreat category above to learn more.",
                    tabMapImageAlt: "Climate and Weather alt.",
                    tabLegendImageAlt: "Legend",
                    tabColor: "var(--color-climate-and-weather)",
                    subThreatPrefix: "CW",
                    subThreatData: [
                        {
                            subThreat: "Change in water temperature",
                            subThreatIcon: "CW-temperature-1979336",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Human-caused climate change.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Temperature affects every aspect of fish physiology, including growth, metabolism, feeding, and reproduction; for example, Atlantic salmon (<em>Salmo salar</em>) show delayed or no ovulation above 13–15°C (<a href='https://documentatiecentrum.watlab.be/imis.php?module=ref&refid=10981' target='_blank'>King and Pankhurst 2000</a>). As waters warm, fish generally shift poleward or upslope to track suitable temperatures (<a href='https://onlinelibrary.wiley.com/doi/10.1111/j.1365-2486.2010.02214.x' target='_blank'>Melles et al. 2011</a>): species with narrow thermal ranges contract, while those with broad ranges expand (<a href='https://doi.org/10.1016/j.biocon.2016.08.021' target='_blank'>Ruiz-Navarro et al. 2016</a>). Warming also changes the spread and virulence of fish diseases and the toxicity of pollutants.",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Industry, fossil fuel sectors, governments, and fishery managers.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Reduce greenhouse gas emissions and increase shade along rivers. Fisheries can adapt by moving fish to locations within their thermal range, increasing disease surveillance, diversifying livelihoods, and taking advantage of new species at higher latitudes.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Fraction of river fish species potentially affected by higher maximum water temperatures under 4.5°C of warming, assuming species cannot move to new areas (<a href='https://doi.org/10.1038/s41467-021-21655-w' target='_blank'>Barbarossa et al. 2021</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Drought",
                            subThreatIcon: "CW-drought-002D5E",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "More frequent and intense periods of abnormally low precipitation.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Less rainfall or snowmelt reduces stream flow and water volume and warms water, degrading water quality and habitat and cutting off connectivity. Fish face more stress, energy deficits, restricted movement, and disease, which can raise mortality, reduce recruitment, and shift community composition, sometimes favoring invasive species (<a href='https://doi.org/10.1007/s11160-018-09545-9' target='_blank'>Lennox et al. 2019</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Industry, fossil fuel sectors, governments, and fishery managers.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Reduce greenhouse gas emissions and water use, remove barriers strategically, and restore habitat. In California, drought compounds existing threats to river trout from habitat change and water withdrawals, making drought-aware water management essential (<a href='https://doi.org/10.1007/s11160-018-09545-9' target='_blank'>Lennox et al. 2019</a>).",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Fraction of river fish species potentially affected by changes in minimum flows under 4.5°C of warming, assuming species cannot move to new areas (<a href='https://doi.org/10.1038/s41467-021-21655-w' target='_blank'>Barbarossa et al. 2021</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Change in flooding",
                            subThreatIcon: "CW-flooding-1760091",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Changes in the frequency and intensity of precipitation events, and sea-level rise.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Flood-borne sediment can deplete oxygen, affect water quality, or destroy nests and eggs. Rising sea-levels affect freshwater salinization. Floods can also breach barriers that contain non-native fish. For example, flooding in the Illinois River could release invasive carp into the Laurentian Great Lakes, threatening a $7 billion fishing industry (<a href='https://www.ideals.illinois.edu/items/92341/bitstreams/301346/data.pdf' target='_blank'>Spacapan et al. 2016</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Water resource and fishery management sectors.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Identify and maintain complex aquatic habitats that give fish refuge during and after floods.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Fraction of river fish species potentially affected by changes in maximum flows under 4.5°C of warming, assuming species cannot move to new areas (<a href='https://doi.org/10.1038/s41467-021-21655-w' target='_blank'>Barbarossa et al. 2021</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Change in wind patterns",
                            subThreatIcon: "CW-wind-7479181",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Warmer water may reduce land–water temperature contrasts, warming oceans and atmosphere may make daily winds more variable, and melting Arctic ice may alter the jet stream, increasing the variability or severity of winds.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Wind is a key driver of when lakes mix after summer stratification into warm surface and cold bottom layers. Earlier and longer stratification can deplete oxygen in productive lakes, potentially increasing fish die-offs.",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Change in wind patterns sectors",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Reduce greenhouse gas emissions and external nutrient loading to maintain oxygen or consider artificial circulation.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Projected change in the number of days lakes are stratified between 2006–2020 and 2041–2060, averaged across four climate models (<a href='https://doi.org/10.48364/ISIMIP.931371.1' target='_blank'>ISIMIP2b Global Lakes Sector</a>; <a href='https://doi.org/10.5281/zenodo.4617643' target='_blank'>Woolway 2021</a>). Data are available for Northern Hemisphere lakes only.",
                                    }],activeOnLoad: false
                                }
                            ]
                        },
                        {
                            subThreat: "Change in ice cover",
                            subThreatIcon: "CW-ice-cover-5314745",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Rising global temperatures.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Less ice cover exposes overwintering eggs to more turbulence and storms, increasing mortality, particularly in poor spawning habitat. For example, ice cover is an important predictor of lake whitefish (<em>Coregonus clupeaformis</em>) recruitment, which is projected to decline as ice cover is lost (<a href='https://doi.org/10.1016/j.jglr.2015.03.015' target='_blank'>Lynch et al. 2015</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Industry and fishery managers.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Reduce greenhouse gas emissions and designate spawning refuges with optimal habitat.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Projected change in the number of days lakes are ice-covered between 2001–2020 and 2041–2060, weighted by lake area and averaged across four climate models (<a href='https://doi.org/10.6084/m9.figshare.19424801.v1' target='_blank'>Wang and Feng 2022</a>). Data are available for Northern Hemisphere lakes only.",
                                    }],activeOnLoad: false
                                }
                            ]
                        }
                    ]
                },
                {
                    tabTitle: "Invasive species",
                    tabContentTitle: "Invasive species",
                    tabContentTitleID: "invasive-species",
                    tabIcon: "invasive-species-973991",
                    tabText: "Invasive non-native species reach new waters accidentally or deliberately. Once established, they compete with and prey on native fish, spread disease, and disrupt food webs. Select a subthreat category above to learn more.",
                    tabMapImageAlt: "Invasive species alt.",
                    tabLegendImageAlt: "Legend",
                    tabColor: "var(--color-invasive-species)",
                    subThreatPrefix: "IS",
                    subThreatData: [
                        {
                            subThreat: "Invasive non-native species",
                            subThreatIcon: "invasive-species-973991",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Accidental introductions occur through aquatic plant transfers, recreational boats, equipment, and water moved between waters, and escapes from aquarium-trade hatcheries during floods or drainage. Invasive fingerlings can also be mixed into fish stocked across borders. Deliberate introductions for aquaculture, stocking of closed water bodies, and biological control can lead to escapes, and non-native fish are sometimes released into open waters during ceremonies or when they outgrow aquariums.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "Once established, invasive species compete with and prey on native species, spread diseases to which native fish have no resistance, and act as parasites or egg predators. Invasive plants such as water hyacinth (<em>Pontederia crassipes</em>) clog waterways, and filter-feeding invaders form dense beds that deplete plankton and can promote weed growth or toxic algal blooms. Together these changes alter nutrient cycling and food webs and can harm commercial fisheries.",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Hatchery operations, particularly smaller hatcheries trading across borders.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Risk assessment and stronger biosecurity practices, particularly when recreational boats are moved between waters.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Area-weighted number of harmful freshwater invasive species in each basin, based on counts by freshwater ecoregion (<a href='https://databasin.org/datasets/11a22d0ce62142258baaa1cbaf520fb2/' target='_blank'>Boucher 2009</a>; <a href='https://esajournals.onlinelibrary.wiley.com/doi/10.1890/070064' target='_blank'>Molnar et al. 2008</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        }
                    ]
                },
                {
                    tabTitle: "Fishing pressure",
                    tabContentTitle: "Fishing pressure",
                    tabContentTitleID: "fishing-pressure",
                    tabIcon: "fishing-pressure-4847095",
                    tabText: "Overfishing occurs when fish are harvested faster than populations can replenish. Although often seen as the main threat to fisheries, fishing pressure ranks lowest of the five threat categories globally. Select a subthreat category above to learn more.",
                    tabMapImageAlt: "Fishing pressure alt.",
                    tabLegendImageAlt: "Legend",
                    tabColor: "var(--color-fishing-pressure)",
                    subThreatPrefix: "FP",
                    subThreatData: [
                        {
                            subThreat: "Overfishing",
                            subThreatIcon: "fishing-pressure-4847095",
                            subThreatText: [
                                // {
                                    // heading: "Drivers",
                                    // content: [{type: "text", content: "Rising demand for fish from large rural populations, local food insecurity, and commercial markets. Many inland fisheries are open access, with little regulation or control of fishing effort, and intense fishing pressure and illegal gear are common. Recreational fisheries can also have high exploitation rates and target large species.",
                                    // }],activeOnLoad: true
                                // },
                                {
                                    heading: "Potential impacts",
                                    content: [{type: "text", content: "As target stocks decline, fishers shift gears and species (<a href='https://doi.org/10.1641/0006-3568(2005)055[1041:OOIW]2.0.CO;2' target='_blank'>Allan et al. 2005</a>). Large, slow-growing species are depleted first and replaced by smaller species lower in the food web (<a href='https://doi.org/10.1006/jmsc.2000.0738' target='_blank'>Murawski 2000</a>). In Lake Victoria, juvenile fish dominated catches in the early 2000s, and catch per boat fell from 80 to 45 kilograms per day (<a href='https://doi.org/10.1111/j.1365-2400.2006.00477.x' target='_blank'>Matsuishi et al. 2006</a>).",
                                    }],activeOnLoad: true
                                },
                                // {
                                    // heading: "Sectors involved",
                                    // content: [{type: "text", content: "Local and national fisheries management bodies, fishers, and the agricultural sector.",
                                    // }],activeOnLoad: false
                                // },
                                {
                                    heading: "Mitigation and adaptation opportunities",
                                    content: [{type: "text", content: "Adopt an ecosystem approach to fisheries management that conserves ecosystem function and trophic structure; introduce regulations such as closed seasons, gear restrictions, and size limits; and support alternative livelihoods, stocking, or aquaculture.",
                                    }],activeOnLoad: false
                                },
                                {
                                    heading: "How the threat is measured",
                                    content: [{type: "text", content: "Documented fish catch in major river basins and lakes important to inland fisheries, allocated to basins by area, river discharge, and human population (<a href='https://doi.org/10.4060/cb2827en' target='_blank'>FAO 2021</a>).",
                                    }],activeOnLoad: false
                                }
                            ]
                        }
                    ]
                }
            ],
        }
    }
}