<template>
    <header id="NavBar" class="header-nav"  role="banner">
        <div class="topnav">
            <RouterLink :to="{ path: '/' }">Home</RouterLink>
            <div class="dropdown">
            <a href="javascript:void(0)" class="dropdown">Projects</a>
            <div class="dropdown-content">
                <RouterLink :to="{ path: '/beaufort-sea' }">Beaufort Sea</RouterLink>
                <RouterLink :to="{ path: '/fire-in-ice' }">Fire in Ice</RouterLink>
                <RouterLink :to="{ path: '/findex' }">Findex</RouterLink>
            </div>
        </div>
            <a href="https://labs.waterdata.usgs.gov/visualizations/index.html#/" target="_blank">About VizLab</a>
        </div>
    </header>
</template>

<script setup>
</script>

<style scoped>

.header-nav {
    position: -webkit-sticky; /* For Safari */
    position: sticky;
    top: 0;
    z-index: 1000;
}

.topnav {
    background-color: #6b7c94;
    overflow: hidden;
    display: flex;
    justify-content: center;
}

/* Style the links inside the navigation bar */
.topnav a {
    color: #f2f2f2;
    text-align: center;
    padding: 14px 16px;
    text-decoration: none;
    font-size: 17px;

}

/* Change the color of links on hover */
.topnav a:hover {
    background-color: #d2e9ff;
    color: var(--dodger-blue);
}

.dropdown {
    float: left;
    overflow: hidden;
}

/* Dropdown button */
.dropdown .dropbtn {
    font-size: 16px;  
    border: none;
    outline: none;
    color: white;
    padding: 14px 16px;
    background-color: inherit;
    font-family: inherit;
    margin: 0;
}

/* Add a background color to dropdown links on hover */
.dropdown-content a:hover {
    background-color: #d2e9ff;
    color: var(--dodger-blue)
}

/* Show the dropdown menu on hover */
.dropdown:hover .dropdown-content {
    display: block;
}

/* Dropdown content (hidden by default) */
.dropdown-content {
    display: none;
    position: absolute;
    top: 48px;
    background-color: #f9f9f9;
    min-width: 160px;
    box-shadow: 0px 8px 16px 0px rgba(0,0,0,0.2);
    z-index: 1;
}

/* Links inside the dropdown */
.dropdown-content a {
    float: none;
    color: black;
    padding: 12px 16px;
    text-decoration: none;
    display: block;
    text-align: left;
}
</style>
