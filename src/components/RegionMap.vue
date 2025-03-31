<template>
    <div ref="mapContainer" class="map-container"></div>
</template>
  
<script setup>
  import { onMounted, ref, watch } from 'vue'
  import * as d3 from 'd3'
  import * as topojson from 'topojson-client'

  const publicPath = import.meta.env.BASE_URL; // this gets the base url for the site

  // S3 resource sourcing
  const s3ProdURL = import.meta.env.VITE_APP_S3_PROD_URL;
  
  const mapContainer = ref(null)
  let mapLayers;

  const emit = defineEmits(['regionSelected']); 

  // props definition, allowing customized paths and datasets
  const props = defineProps({
    layerVisibility: {
      type: Object,
      required: true
    },
    layerPaths: {
      type: Object,
      required: true
    },
    layerMag: {
      type: String,
      required: true
    },
    layerX: {
      type: String,
      required: true
    },
    layerY: {
      type: String,
      required: true
    },
    regionsDataUrl: {
      type: String,
      required: true
    },
    regionsVar: {
      type: String,
      required: true
    },
    usOutlineUrl: {
      type: String,
      required: true
    }
  })

  const updateLayers = () => {
    if (!mapLayers) return 

    const visibleLayers = Object.entries(props.layerVisibility).map(([key, value]) => ({
      key,
      path: props.layerPaths[key]?.path,
      visible: value
    }))

    mapLayers.selectAll('image')
      .data(visibleLayers, d => d.key) // use key as the identifier
      .join(
        enter => enter.append('image')
          .attr('xlink:href', d => s3ProdURL + 'images/water-availability/' + d.path)
          .attr('x', props.layerX)
          .attr('y', props.layerY)
          .attr('width', 800 * props.layerMag)
          .attr('height', 550 * props.layerMag)
          .style('opacity', d => (d.visible ? 1 : 0))
          .style('display', d => (d.visible ? 'block' : 'none')), // ensure full hiding
        update => update
          .style('opacity', d => (d.visible ? 1 : 0))
          .style('display', d => (d.visible ? 'block' : 'none')), // ensure full hiding
        exit => exit.remove()
      )
      
  }

  // watch layerVisibility for changes
  watch(
    () => props.layerVisibility,
    () => {
      updateLayers(); // Trigger layer updates
    },
    { deep: true }
  );

  
  onMounted(async () => {
    if (!mapContainer.value) {
      console.error('mapContainer is not defined')
      return
    }
  
    const width = 800;
    const height = 550;
    const maxHeight = 800;
  
    // create svg that holds the map
    const svg = d3.select(mapContainer.value)
      .append('svg')
      .attr('viewBox', `0 0 ${width} ${height}`)
      .attr('preserveAspectRatio', 'xMidYMid meet')
      .classed('responsive-svg', true);

    mapLayers = svg.append('g').attr('class', 'map-layers')
    updateLayers()  
  
    // resizing so flexes to page width but stays within reasonable height
    const resizeSvg = () => {
      const containerWidth = mapContainer.value.clientWidth;
      const containerHeight = Math.min(mapContainer.value.clientHeight, maxHeight);
      svg.attr('width', containerWidth).attr('height', containerHeight);
    };
  
    resizeSvg();
    window.addEventListener('resize', resizeSvg);
  
    try {
    // read in data
      // region shapes - feature collection
      const topoRegions = await d3.json(`${publicPath}${props.regionsDataUrl}`);
      const geoRegions = topojson.feature(topoRegions, topoRegions.objects[Object.keys(topoRegions.objects)[0]]);

      // CONUS outline - single feature
      const topoUS = await d3.json(props.usOutlineUrl);
      const geoUS = topojson.feature(topoUS, topoUS.objects['foo']);  
      const projection = d3.geoIdentity().reflectY(true).fitSize([width, height], geoRegions);
      const path = d3.geoPath().projection(projection);
  
      // Overlay the raster images
      const scale_size = 1.2; // scaling pngs because they have an added margin when exported from ggplot
      svg.append('g')
        .selectAll('image')
        .data(props.layerPaths)
        .enter()
        .append('image')
        .attr('xlink:href', d => import.meta.env.BASE_URL + d)
        .attr('x', -80) // nudging png to fit within svg bounds
        .attr('y', -55) // nudging png to fit within svg bounds
        .attr('width', width * scale_size)
        .attr('height', height * scale_size);
  

      // draw region boundaries
      svg.append('g')
        .selectAll('path')
        .data(geoRegions.features)
        .join('path')
        .attr('d', path)
        .attr('class', d => `region ${d[props.RegionsVar]}`)
        .attr('fill', 'transparent')
        .attr("opacity", 0)
        .attr('stroke', 'black')
        .attr('stroke-width', '1px')
        // add interaction to highlight selected region and update bar chart
        .on('mouseover',function(event, d) {
            d3.selectAll('.region')
                .attr('fill','lightgrey')
                .attr('opacity', 0.8)
                .attr('stroke',"black")

            // highlight the selected region with transparent fill
            d3.select(this)
                .attr('fill', 'transparent')
                .attr('opacity', 1)
                .attr('stroke', 'black')
                .attr('stroke-width', '1.5px')
                .raise(); // bring the selected region to the front

            // update the bar chart with the selected region's data
             highlightRegionAndUpdateChart(event, d);

        })
        .on('mouseout', function () {
            // reset all regions to transparent fill and opacity 0
            d3.selectAll('.region')
            .attr('fill', 'transparent')
            .attr('opacity', 0)
            .attr('stroke', 'white')
            .attr('stroke-width', '1.2px');

            // reset bar chart to default aggregated data
            emit('regionSelected', 'lower 48 United States');
        })
        .on('click', (event, d) => {
          d3.selectAll('.region')
                .attr('fill','lightgrey')
                .attr('opacity', 0.8)
                .attr('stroke',"black")

            // highlight the selected region with transparent fill
            d3.select(this)
                .attr('fill', 'transparent')
                .attr('opacity', 1)
                .attr('stroke', 'white')
                .attr('stroke-width', '1.5px')
                .raise(); // bring the selected region to the front
          highlightRegionAndUpdateChart(event, d);
        });

    // add double outline for CONUS
    svg.append('g')
        .append('path')
        .datum(geoUS)
        .attr("class", "outline-conus")
        .attr('d', path)
        .attr('fill', 'none')
        .attr('stroke', 'grey')
        .attr('stroke-width', '2px');

      // add another outline to regions to make it more visible? idk if this looks great
      svg.append('g')
        .selectAll('path')
        .data(geoRegions.features)
        .join('path')
        .attr('d', path)
        .attr('fill', 'none')
        .attr('stroke', 'white')
        .attr('stroke-width', '0.75px');


        d3.selectAll('.region')
                .attr('fill','lightgrey')
                .attr('opacity', 0.8)
                .attr('stroke',"black")

        // highlight the selected region with transparent fill
        d3.select(this)
            .attr('fill', 'transparent')
            .attr('opacity', 1)
            .attr('stroke', 'black')
            .attr('stroke-width', '1.5px')
  
    } catch (error) {
      console.error('Error loading TopoJSON:', error);
    }
  });

  // selection effects and filtering with interaction
  function highlightRegionAndUpdateChart(event, d) {
    // update bar chart with regional data
    const regionClassFilter = d.properties.Region_nam;
    emit('regionSelected', regionClassFilter); // send region_nam to parent
  }
</script>
  
<style>
  .map-container {
    display: flex;
    justify-content: center;
    align-items: center;
    width: 100%;
    height: auto;
    max-height: 700px;
  }
  
  .responsive-svg {
    width: 100%;
    height: auto;
    max-height: 100%;
  }
  
  .bar-chart-svg {
    width: 100%;
    height: auto;
    max-height: 100%;
  }
  .outline-conus {
    filter: drop-shadow(0px 0px 0px rgba(2, 2, 2, 0.5));
  }
</style>
  