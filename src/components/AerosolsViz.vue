<template>
    <!---VizSection-->
    <VizSection
        id="cross-section"
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <template #heading>
            <h2>
                {{ text.heading }}
            </h2>
        </template>
        <!-- FIGURES -->
        <template #aboveExplanation>
            <p v-html="text.paragraph1" />
        </template>
        <template #figures>
            <div id="aerosols-grid-container" />
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
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {
        try {
            // Use external svg from s3
            d3.xml("https://labs.waterdata.usgs.gov/visualizations/svgs/aerosols_map.svg").then(function(xml) {
                // add svg content to DOM
                const svgGrid = document.getElementById("aerosols-grid-container")
                svgGrid.appendChild(xml.documentElement);
                
                // add id to svg
                d3.select("#aerosols-grid-container").select("svg")
                    .attr("id", "aerosols-svg")

                // add interactivity
                addInteractions();
            });
        } catch (error) {
            console.error('Error during component mounting', error);
        }        
    });


    function draw_trajectories(fire,trajectory_num,smoke_lines,delay,duration){
        for(let smoke_num=0;smoke_num<smoke_lines;smoke_num++){
            var trajectory_line = d3.select("#trajectory-" + fire + "-" + trajectory_num + "-" + smoke_num).selectAll("path");
            var totalLength = trajectory_line.node().getTotalLength();
            var opacity = 0.6;
            
            // Animate the line
            trajectory_line.style("transition", "none")
                .style("stroke-dasharray", totalLength + " " + totalLength)
                .style("stroke-dashoffset", totalLength - 1) //minus 1 remove artifacts that show before animation.
                .style("stroke-opacity", opacity)
                .transition()
                .delay(delay)
                .duration(duration)
                .ease(d3.easeLinear)
                .style("stroke-dashoffset", 0);
        }
    }

    function remove_trajectories(fire,trajectory_num,smoke_lines){
        for(let smoke_num=0;smoke_num<smoke_lines;smoke_num++){
            d3.select("#trajectory-" + fire + "-" + trajectory_num + "-" + smoke_num).selectAll("path")
                .style("stroke-opacity", 0.0)
                .interrupt();
        }
    }

    function mouseover(event,number_of_fires,smoke_lines) {
        if (event.currentTarget.id.startsWith("source-")){
            for (let fire=0;fire<number_of_fires;fire++){
                    d3.select('#wildfire-label-'+fire).selectAll("path")
                        .style("opacity", 0.0);
                    d3.select('#wildfire-label-'+fire).selectAll("text")
                        .style("opacity", 0.0);
            }
            d3.select('#multi-path-label').selectAll("text")
                .style("opacity", 0.0);
            d3.select('#multi-path-label').selectAll("path")
                .style("opacity", 0.0);
            let fire = event.currentTarget.id.slice(7);
            for(let trajectory_num=1;trajectory_num<=24;trajectory_num++){
                draw_trajectories(fire,trajectory_num,smoke_lines,0.0,2700);//300*(trajectory_num-1),2700);
            }
        }
    }

    function mouseout(event,number_of_fires,smoke_lines) {
        if (event.currentTarget.id.startsWith("source-")){
            for (let fire=0;fire<number_of_fires;fire++){
                    d3.select('#wildfire-label-'+fire).selectAll("path")
                        .style("opacity", 0.75);
                    d3.select('#wildfire-label-'+fire).selectAll("text")
                        .style("opacity", 1.0);
            }
            d3.select('#multi-path-label').selectAll("path")
                .style("opacity", 0.75);
            d3.select('#multi-path-label').selectAll("text")
                .style("opacity", 1.0);
            let fire = event.currentTarget.id.slice(7);
            for(let trajectory_num=1;trajectory_num<=24;trajectory_num++){
                remove_trajectories(fire,trajectory_num,smoke_lines);
            }
        }
    }

    function mouseenter(event,default_fire,number_of_fires,default_smoke_num,smoke_lines) {
        if (event.currentTarget.id.startsWith("figure_1")){
            d3.select('#single-path-label').selectAll("path")
                .style("opacity", 0.0);
            d3.select('#single-path-label').selectAll("text")
                .style("opacity", 0.0);
            d3.select('#multi-path-label').selectAll("path")
                .style("opacity", 0.75);
            d3.select('#multi-path-label').selectAll("text")
                .style("opacity", 1.0);
            remove_trajectories(default_fire,default_smoke_num,smoke_lines)
        }
    }

    function mouseleave(event,default_fire,number_of_fires,default_smoke_num,smoke_lines) {
        if (event.currentTarget.id.startsWith("figure_1")){
            d3.select('#single-path-label').selectAll("path")
                .style("opacity", 0.75);
            d3.select('#single-path-label').selectAll("text")
                .style("opacity", 1.0);
            d3.select('#multi-path-label').selectAll("text")
                .style("opacity", 0.0);
            d3.select('#multi-path-label').selectAll("path")
                .style("opacity", 0.0);
            draw_trajectories(default_fire,default_smoke_num,smoke_lines,0,1000)
        }
    }

    function addInteractions() {
        // set viewbox for svg with loss function chart
        const aerosolsSVG = d3.select("#aerosols-svg");

        // draw default line
        const default_fire = 1;
        const default_smoke_num = 8;
        const smoke_lines = 1;
        const number_of_fires = 2;
        for (let fire=0;fire<number_of_fires;fire++){
            d3.select('#wildfire-label-'+fire).selectAll("path")
                .style("opacity", 0.75);
            d3.select('#wildfire-label-'+fire).selectAll("text")
                .style("opacity", 1.0);
        }
        d3.select('#JIF-label').selectAll("path")
            .style("opacity", 0.75);
        d3.select('#JIF-label').selectAll("text")
            .style("opacity", 1.0);
        d3.select('#single-path-label').selectAll("path")
            .style("opacity", 0.75);
        d3.select('#single-path-label').selectAll("text")
            .style("opacity", 1.0);
        d3.select('#multi-path-label').selectAll("path")
            .style("opacity", 0.0);
        d3.select('#multi-path-label').selectAll("text")
            .style("opacity", 0.0);

        draw_trajectories(default_fire,default_smoke_num,smoke_lines,0,1000)

        // Add interaction to loss function chart
        aerosolsSVG.selectAll("g")
            .on("mouseover", (event) => mouseover(event,number_of_fires,smoke_lines))
            .on("mouseout", (event) => mouseout(event,number_of_fires,smoke_lines))
            .on("mouseenter", (event) => mouseenter(event,default_fire,number_of_fires,default_smoke_num,smoke_lines))
            .on("mouseleave", (event) => mouseleave(event,default_fire,number_of_fires,default_smoke_num,smoke_lines));
    }
</script>

<style scoped lang="scss">
    #aerosols-grid-container {
        display: grid;
        width: 100%;
        max-width: 1200px;
        margin: 0 auto 0 auto;
        grid-template-areas:
            "chart";
    }
</style>

<style lang="scss">
/* css for elements added/classed w/ d3 */
    #aerosols-svg {
        grid-area: chart;
        place-self: center;
        max-height: 100%;
        max-width: 100%;
    }
</style>

<style>
    #wildfire-label-0 {
        cursor: default;
    }
    #wildfire-label-1 {
        cursor: default;
    }
    #single-path-label {
        cursor: default;
    }
    #multi-path-label {
        cursor: default;
    }
</style>