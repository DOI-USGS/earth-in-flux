<template>
    <div class="chart">
        <div class="takeaway">
            <p>{{ description }}</p>
        </div>
        <img :src="src" :alt="alt"/>        
        <!--div class="project">
            <p>{{ project }}</p>
        </div-->
    </div>
</template>

<script setup>
    defineProps({
        src: {
            type: String,
            default: ``
        },
        alt: {
            type: String,
            default: ``
        },
        project: {
            type: String,
            default: ''
        },
        description: {
            type: String,
            default: ''
        }
    })
</script>

<style scoped lang="scss">
    $collapsed-height: 280px;
    $expanded-height: 400px;
    $height-difference: $expanded-height - $collapsed-height;
	.chart{
		display: flex;
		flex-direction: column-reverse;
		justify-content: flex-start;
		position: relative;
		width: $collapsed-height;
		height: $collapsed-height;
        -moz-transition: height 1s ease;
        -webkit-transition: height 1s ease;
        -o-transition: height 1s ease;
        transition: height 1s ease;
        background-color: var(--faded-usgs-blue);
        box-shadow: 0px 0px 8px rgba(39,44,49,.07), 1px 4px 4px rgba(39,44,49,.04);
	}
    .chart:hover {
        height: $expanded-height;
        transform: scaleY(1);
    }
    .expanded {
        height: $expanded-height;
        transform: scaleY(1);
    }
    .takeaway {
        position: absolute;
        left: 0;
        right: 0;
        bottom: calc($height-difference / 2);
        line-height: 2rem;
        -webkit-transition: top 1s ease; 
        -moz-transition: top 1s ease; 
        -ms-transition: top 1s ease; 
        -o-transition: top 1s ease; 
        transition: top 1s ease;
    }
    .takeaway p {
        font-family: sans-serif; /* This is fallback font for old browsers */
        font-family: var(--title-font);
        color: var(--color-title-text);
        transform: translateY(50%);
        padding: 0.2rem 2rem 0.2rem 2rem;
        opacity: 0;
        transition: all 1s ease; 
    }
    .chart:hover .takeaway p {
        opacity: 1;
    }
    .expanded .takeaway p {
        opacity: 1;
    }
    .chart img {
        position: absolute;
        width: 100%;
        top: 0;
        background-color: white;
        transition: all 1s ease; 
    }
    .chart:hover img {
        transform: translate3D(0,-0.5px,0) scale(1.05); 
        box-shadow:  0px 0px 8px rgba(39,44,49,.07), 4px 8px 4px rgba(39,44,49,.1), 1px 4px 4px rgba(39,44,49,.07);
    }
    .project {
        font-family: sans-serif; /* This is fallback font for old browsers */
        font-family: var(--title-font);
        color: var(--color-title-text);
        position: absolute;
        top: $collapsed-height;
        width: auto;
        left: 50%;
        transform: translate(-50%, -50%);
        -webkit-transition: top 1s ease; 
        -moz-transition: top 1s ease; 
        -ms-transition: top 1s ease; 
        -o-transition: top 1s ease; 
        transition: top 1s ease;
        background-color: var(--faded-usgs-blue);
        border-radius: 0.25rem 0.25rem 0 0;
        box-shadow: 1px 4px 4px rgba(107, 107, 107, 0.1);
    }
    .project p {
        font-size: 1.4rem;
        padding: 0.9rem 1.5rem 0.9rem 1.5rem; 
        text-align: center;
        width: auto;
        font-weight: 400;
    }
    .chart:hover .project {
        top: $expanded-height;
    }
    .expanded .project {
        top: $expanded-height;
    }
</style>