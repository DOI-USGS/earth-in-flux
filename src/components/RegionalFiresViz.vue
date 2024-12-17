<template>
    <!---VizSection-->
    <VizSection
        id="cross-section"
        :figures="true"
        :fig-caption="false"
    >
        <template #heading>
            <h2 />
        </template>
        <!-- FIGURES -->
        <template #aboveExplanation>
            <p v-html="text.paragraph1" />
            <p v-html="text.paragraph2" />
        </template>
        <template #figures>
            <div id="aerosols-grid-container" />
            <div class="text-container">
                <button id="reset-button" @click="resetViz">Reset map</button>
            </div>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
        </template>
    </VizSection>
</template>

<script setup>
    import { onMounted } from "vue";
    import * as d3 from 'd3';
    import { isMobile } from 'mobile-device-detect';
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })

    // global variables
    const mobileView = isMobile;
    const default_fire = 1;
    const default_smoke_num = 8;
    const smoke_lines = 1;
    const number_of_fires = 2;
    const smoke_opacity = 0.6;
    const smoke_width = 15;

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {
        try {
            // Use external svg from s3
            d3.xml("https://labs.waterdata.usgs.gov/visualizations/svgs/regional_fires_map_v10.svg").then(function(xml) {
                // add svg content to DOM
                const svgGrid = document.getElementById("aerosols-grid-container")
                svgGrid.appendChild(xml.documentElement);
                
                // add id to svg
                const aerosolsSVG = d3.select("#aerosols-grid-container").select("svg")
                    .attr("id", "aerosols-svg")

                aerosolsSVG
                    .append("title")
                    .text("A shaded-relief map of southeastern Alaska, showing the location of the Juneau Icefield. The map also shows the location of the 2015 Dennison Fork wildfire and the 2016 Steamboat Creek wildfire, both northwest of the icefield. Transparent grey lines branching out from each wildfire location show modeled paths of smoke transport from these fires. Some of the potential smoke paths pass over the Juneau Icefield, suggesting that aerosols from these fires may have been deposited on the icefield.");

                // add interactivity
                addInteractions();
            });
        } catch (error) {
            console.error('Error during component mounting', error);
        }        
    });


    function draw_trajectories(fire,trajectory_num,delay,duration){
        for(let smoke_num=0;smoke_num<smoke_lines;smoke_num++){
            const trajectory_line = d3.select("#trajectory-" + fire + "-" + trajectory_num + "-" + smoke_num).selectAll("path");
            const totalLength = trajectory_line.node().getTotalLength();
            
            // Animate the line
            trajectory_line.style("transition", "none")
                .style("filter", "url(#glow)")
                .style("stroke-dasharray", totalLength + " " + totalLength)
                .attr("stroke-dashoffset", totalLength - 1) //minus 1 remove artifacts that show before animation.
                .style("stroke-opacity", smoke_opacity)
                .style('stroke-linecap', 'round')
                .style("stroke-width", smoke_width)
                .transition()
                .delay(delay)
                .duration(duration)
                .ease(d3.easeLinear)
                .attr("stroke-dashoffset", 0);
        }
    }

    function remove_trajectories(fire,trajectory_num){
        for(let smoke_num=0;smoke_num<smoke_lines;smoke_num++){
            d3.select("#trajectory-" + fire + "-" + trajectory_num + "-" + smoke_num).selectAll("path")
                .style("stroke-opacity", 0.0)
                .interrupt();
        }
    }

    function mouseover(event) {
        if (event.currentTarget.id.startsWith("source-")){
            for (let fire=0;fire<number_of_fires;fire++){
                    d3.select('#wildfire-label-'+fire).selectAll("path")
                        .style("opacity", 0.0);
                    d3.select('#wildfire-label-'+fire).selectAll("text")
                        .style("opacity", 0.0);
            }
            if (mobileView == true){
                d3.select('#multi-path-label-mb-2').selectAll("path")
                    .style("opacity", 0.0);
                d3.select('#multi-path-label-mb-2').selectAll("text")
                    .style("opacity", 0.0);
            }

            const fire = event.currentTarget.id.slice(7);
            let otherFire = fire === '1' ? '0' : '1';
            for(let trajectory_num=1;trajectory_num<=24;trajectory_num++){
                remove_trajectories(otherFire,trajectory_num);
            }
            for(let trajectory_num=1;trajectory_num<=24;trajectory_num++){
                draw_trajectories(fire,trajectory_num,0.0,2700);//300*(trajectory_num-1),2700);
            }
        }
    }

    function mouseout(event) {
        if (event.currentTarget.id.startsWith("source-")){
            for (let fire=0;fire<number_of_fires;fire++){
                    d3.select('#wildfire-label-'+fire).selectAll("path")
                        .style("opacity", 0.75);
                    d3.select('#wildfire-label-'+fire).selectAll("text")
                        .style("opacity", 1.0);
            }

            const fire = event.currentTarget.id.slice(7);
            for(let trajectory_num=1;trajectory_num<=24;trajectory_num++){
                remove_trajectories(fire,trajectory_num);
            }
        }
    }

    function mouseenter(event) {
        if (event.currentTarget.id.startsWith("figure_1")){
            for (let fire=0;fire<number_of_fires;fire++){
                d3.select('#wildfire-label-'+fire).selectAll("path")
                    .style("opacity", 0.75);
                d3.select('#wildfire-label-'+fire).selectAll("text")
                    .style("opacity", 1.0);
                for(let trajectory_num=1;trajectory_num<=24;trajectory_num++){
                    remove_trajectories(fire, trajectory_num);
                }
            }
            d3.select('#single-path-label').selectAll("path")
                .style("opacity", 0.0);
            d3.select('#single-path-label').selectAll("text")
                .style("opacity", 0.0);
            if (mobileView == true){
                d3.select('#multi-path-label-mb-1').selectAll("path")
                    .style("opacity", 0.0);
                d3.select('#multi-path-label-mb-1').selectAll("text")
                    .style("opacity", 0.0);
                d3.select('#multi-path-label-mb-2').selectAll("path")
                    .style("opacity", 0.75);
                d3.select('#multi-path-label-mb-2').selectAll("text")
                    .style("opacity", 1.0);
            } else{
                d3.select('#multi-path-label-dt').selectAll("path")
                    .style("opacity", 0.0);
                d3.select('#multi-path-label-dt').selectAll("text")
                    .style("opacity", 0.0);
            }
            remove_trajectories(default_fire,default_smoke_num)
        }
    }

    function mouseleave(event) {
        if (event.currentTarget.id.startsWith("figure_1")){
            for (let fire=0;fire<number_of_fires;fire++){
                d3.select('#wildfire-label-'+fire).selectAll("path")
                    .style("opacity", 0.0);
                d3.select('#wildfire-label-'+fire).selectAll("text")
                    .style("opacity", 0.0);
            }
            d3.select('#single-path-label').selectAll("path")
                .style("opacity", 0.75);
            d3.select('#single-path-label').selectAll("text")
                .style("opacity", 1.0);
            if (mobileView == true){
                d3.select('#multi-path-label-mb-1').selectAll("path")
                    .style("opacity", 0.75);
                d3.select('#multi-path-label-mb-1').selectAll("text")
                    .style("opacity", 1.0);
                d3.select('#multi-path-label-mb-2').selectAll("path")
                    .style("opacity", 0.0);
                d3.select('#multi-path-label-mb-2').selectAll("text")
                    .style("opacity", 0.0);
            } else{
                d3.select('#multi-path-label-dt').selectAll("path")
                    .style("opacity", 0.75);
                d3.select('#multi-path-label-dt').selectAll("text")
                    .style("opacity", 1.0);
            }
            draw_trajectories(default_fire,default_smoke_num,0,1000)
        }
    }

    function addInteractions() {
        // set width and height of svg
        const aerosolsSVG = d3.select("#aerosols-svg")
            .attr("width", "100%")
            .attr("height", "100%");

        // Container for the gradients
        const defs = aerosolsSVG.append("defs");

        // Filter for the outside glow
        const filter = defs.append("filter")
            .attr("id","glow");
        // wide blur
        filter.append("feGaussianBlur")
            .attr("stdDeviation","1.5")
            .attr("result","coloredBlur");

        // Set display of all labels, and hide from screen reader
        for (let fire=0;fire<number_of_fires;fire++){
            d3.select('#wildfire-label-'+fire).selectAll("path")
                .style("opacity", 0.0);
            d3.select('#wildfire-label-'+fire).selectAll("text")
                .attr("class", "chart-text")
                .style("opacity", 0.0)
                .attr("aria-hidden", true);
        }
        d3.select('#JIF-label').selectAll("path")
            .style("opacity", 0.75);
        d3.select('#JIF-label').selectAll("text")
            .attr("class", "chart-text")
            .style("opacity", 1.0)
            .attr("aria-hidden", true);
        d3.select('#single-path-label').selectAll("path")
            .style("opacity", 0.75);
        d3.select('#single-path-label').selectAll("text")
            .attr("class", "chart-text")
            .style("opacity", 1.0)
            .attr("aria-hidden", true);
        d3.select('#multi-path-label-mb-2').selectAll("path")
            .style("fill", "#151515")
            .style("fill-opacity", 0.75)
            .style("opacity", 0.0);
        d3.select('#multi-path-label-mb-2').selectAll("text")
            .attr("class", "chart-text")
            .style("opacity", 0.0)
            .attr("aria-hidden", true);
        if (mobileView == true){
            d3.select('#multi-path-label-mb-1').selectAll("path")
                .style("fill", "#151515")
                .style("fill-opacity", 0.75)
                .style("opacity", 0.75);
            d3.select('#multi-path-label-mb-1').selectAll("text")
                .attr("class", "chart-text")
                .style("opacity", 1.0)
                .attr("aria-hidden", true);
            d3.select('#multi-path-label-dt').selectAll("path")
                .style("opacity", 0.0);
            d3.select('#multi-path-label-dt').selectAll("text")
                .attr("class", "chart-text")
                .style("opacity", 0.0)
                .attr("aria-hidden", true);
        } else{
            d3.select('#multi-path-label-dt').selectAll("path")
                .style("opacity", 0.75);
            d3.select('#multi-path-label-dt').selectAll("text")
                .attr("class", "chart-text")
                .style("opacity", 1.0)
                .attr("aria-hidden", true);
            d3.select('#multi-path-label-mb-1').selectAll("path")
                .style("opacity", 0.0);
            d3.select('#multi-path-label-mb-1').selectAll("text")
                .attr("class", "chart-text")
                .style("opacity", 0.0)
                .attr("aria-hidden", true);
        }

        // Hide axis labels from screen reader and assign class for styling
        const idsToHide = [...Array(10).keys()].slice(1, 10);
        idsToHide.forEach(id => {
            const idGroup = aerosolsSVG.select(`#text_${id}`)

            idGroup
                .attr("aria-hidden", true);

            idGroup.selectAll('text').selectAll('tspan')
                .attr("class", "chart-text")

        })

        // add tab-index selection and keypress events for fire markers
        for (let i=0; i < 2; i++) {
            d3.select(`#source-${i}`)
                .attr('tabindex', 0)
                .attr("role", "button")
                .attr("aria-label", i == 0 ? 'Dennison Fork 2015 wildfire' : 'Steamboat Creek 2016 wildfire') 
                .on("keydown", function(event) {
                    if(event.code == 'Enter' | event.code == 'Space'){
                        keypressFire(event)
                    }
                })
        }

        // draw default line
        draw_trajectories(default_fire,default_smoke_num,0,1000)

        // Add interaction
        aerosolsSVG.selectAll("g")
            .on("mouseover", (event) => mouseover(event))
            .on("mouseout", (event) => mouseout(event))
            .on("mouseenter", (event) => mouseenter(event))
            .on("mouseleave", (event) => mouseleave(event));
    }

    function keypressFire(event) {
        const fire = event.currentTarget.id.slice(7);
        let otherFire = fire === '1' ? '0' : '1';
        for(let trajectory_num=1;trajectory_num<=24;trajectory_num++){
            remove_trajectories(otherFire,trajectory_num);
        }
        for (let fire=0;fire<number_of_fires;fire++){
            d3.select('#wildfire-label-'+fire).selectAll("path")
                .style("opacity", 0.75);
            d3.select('#wildfire-label-'+fire).selectAll("text")
                .style("opacity", 1.0);
        }
        d3.select('#single-path-label').selectAll("path")
            .style("opacity", 0.0);
        d3.select('#single-path-label').selectAll("text")
            .style("opacity", 0.0);
        d3.select('#multi-path-label-dt').selectAll("path")
            .style("opacity", 0.0);
        d3.select('#multi-path-label-dt').selectAll("text")
            .style("opacity", 0.0);
        remove_trajectories(default_fire,default_smoke_num)

        for(let trajectory_num=1;trajectory_num<=24;trajectory_num++){
            draw_trajectories(fire,trajectory_num,0.0,2700);
        }
    }

    function resetViz() {
        for (let fire=0;fire<number_of_fires;fire++){
                d3.select('#wildfire-label-'+fire).selectAll("path")
                    .style("opacity", 0.75);
                d3.select('#wildfire-label-'+fire).selectAll("text")
                    .style("opacity", 1.0);
        }

        for(let trajectory_num=1;trajectory_num<=24;trajectory_num++){
            remove_trajectories(0,trajectory_num);
        }
        for(let trajectory_num=1;trajectory_num<=24;trajectory_num++){
            remove_trajectories(1,trajectory_num);
        }
        for (let fire=0;fire<number_of_fires;fire++){
            d3.select('#wildfire-label-'+fire).selectAll("path")
                .style("opacity", 0.0);
            d3.select('#wildfire-label-'+fire).selectAll("text")
                .style("opacity", 0.0);
        }
        d3.select('#single-path-label').selectAll("path")
            .style("opacity", 0.75);
        d3.select('#single-path-label').selectAll("text")
            .style("opacity", 1.0);
        if (mobileView == true){
            d3.select('#multi-path-label-mb-1').selectAll("path")
                .style("opacity", 0.75);
            d3.select('#multi-path-label-mb-1').selectAll("text")
                .style("opacity", 1.0);
            d3.select('#multi-path-label-mb-2').selectAll("path")
                .style("opacity", 0.0);
            d3.select('#multi-path-label-mb-2').selectAll("text")
                .style("opacity", 0.0);
        } else{
            d3.select('#multi-path-label-dt').selectAll("path")
                .style("opacity", 0.75);
            d3.select('#multi-path-label-dt').selectAll("text")
                .style("opacity", 1.0);
        }
        draw_trajectories(default_fire,default_smoke_num,0,1000)
    }
</script>

<style scoped lang="scss">
    #aerosols-grid-container {
        display: grid;
        width: 100%;
        max-width: 1200px;
        max-height: 85vh;
        margin: 0 auto 0 auto;
        grid-template-areas:
            "chart";
    }
    #reset-button {
        margin: 4rem 0 2rem 0;
    }
</style>

<style lang="scss">
/* css for elements added/classed w/ d3 */
    #aerosols-svg {
        grid-area: chart;
        place-self: center;
        height: 100%;
        max-height: calc(min(85vh, 700px));
        max-width: 100%;
    }
</style>

<style>
    #wildfire-label-0 {
        cursor: default;       
        @media only screen and (max-width: 600px) {            
            transform-origin: 50% 30%;
            transform: scale(1.2,1.2);
        }
    }
    #wildfire-label-1 {
        cursor: default; 
        @media only screen and (max-width: 600px) {            
            transform-origin: 10% 50%;
            transform: scale(1.2,1.2);
        }
    }
    #JIF-label {
        cursor: default;
        @media only screen and (max-width: 600px) {            
            transform-origin: 75% 85%;
            transform: scale(1.2,1.2);
        }
    }
    #single-path-label {
        cursor: default;
    }
    #multi-path-label-dt {
        cursor: default;
    }
    #multi-path-label-mb-1 {
        cursor: default;
    }
    #multi-path-label-mb-2 {
        cursor: default;
    }
</style>