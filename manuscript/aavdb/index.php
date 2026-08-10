<?php
#start a secure session
session_set_cookie_params(['secure'=>true,'httponly'=>true,'samesite'=>'Lax']);
session_start();

#hardcoded credentials
$user = 'tap';
$pass = 'slab';
$pass2 = 'reviewer';


#handle login
if ($_SERVER['REQUEST_METHOD']==='POST') {
    if ($_POST['password']===$pass || $_POST['password']===$pass2) {
        $_SESSION['auth']=true;
        header('Location: /');  #refresh to hide POST
        exit;
    } else {
        $error = 'invalid credentials';
    }
}

#handle logout
if (isset($_GET['logout'])) {
    session_destroy();
    header('Location: /');
    exit;
}

#require auth before showing main content
if (!($_SESSION['auth'] ?? false)):
?>
<!doctype html>
<html>
<head>
  <title>Login</title>
  <style>
    body{font-family:sans-serif;display:flex;height:100vh;justify-content:center;align-items:center;background:#f6f6f6;}
    form{background:#fff;padding:2rem;border-radius:1rem;box-shadow:0 2px 8px rgba(0,0,0,0.1);}
    input{display:block;margin:.5rem 0;padding:.5rem;width:200px;}
    button{padding:.5rem 1rem;}
    p{color:red;}
  </style>
</head>
<body>
  <form method="post">
    <h2>Login</h2>
    <?php if(!empty($error)) echo "<p>$error</p>"; ?>
    <input type="password" name="password" placeholder="password" required>
    <button type="submit">Sign in</button>
  </form>
</body>
</html>
<?php
exit; #stop processing after showing login form
endif;
?>

<!doctype html>
<html lang="en" data-bs-theme="auto">
   <head>
      <meta charset="utf-8" />
      <meta name="viewport" content="width=device-width, initial-scale=1" />
      <meta name="description" content="" />
      <title>Tropism Analysis Portal | AAVdb</title>
      <script src="assets/js/color-modes.js"></script>
      <link href="assets/dist/css/bootstrap.min.css" rel="stylesheet" />
      <meta name="theme-color" content="#712cf9" />
      <style>
         .bd-placeholder-img {
         font-size: 1.125rem;
         text-anchor: middle;
         -webkit-user-select: none;
         -moz-user-select: none;
         user-select: none;
         }
         @media (min-width: 768px) {
         .bd-placeholder-img-lg {
         font-size: 3.5rem;
         }
         }
         .b-example-divider {
         width: 100%;
         height: 3rem;
         background-color: #0000001a;
         border: solid rgba(0, 0, 0, 0.15);
         border-width: 1px 0;
         box-shadow:
         inset 0 0.5em 1.5em #0000001a,
         inset 0 0.125em 0.5em #00000026;
         }
         .b-example-vr {
         flex-shrink: 0;
         width: 1.5rem;
         height: 100vh;
         }
         .bi {
         vertical-align: -0.125em;
         fill: currentColor;
         }
         .nav-scroller {
         position: relative;
         z-index: 2;
         height: 2.75rem;
         overflow-y: hidden;
         }
         .nav-scroller .nav {
         display: flex;
         flex-wrap: nowrap;
         padding-bottom: 1rem;
         margin-top: -1px;
         overflow-x: auto;
         text-align: center;
         white-space: nowrap;
         -webkit-overflow-scrolling: touch;
         }
         .btn-bd-primary {
         --bd-violet-bg: #712cf9;
         --bd-violet-rgb: 112.520718, 44.062154, 249.437846;
         --bs-btn-font-weight: 600;
         --bs-btn-color: var(--bs-white);
         --bs-btn-bg: var(--bd-violet-bg);
         --bs-btn-border-color: var(--bd-violet-bg);
         --bs-btn-hover-color: var(--bs-white);
         --bs-btn-hover-bg: #6528e0;
         --bs-btn-hover-border-color: #6528e0;
         --bs-btn-focus-shadow-rgb: var(--bd-violet-rgb);
         --bs-btn-active-color: var(--bs-btn-hover-color);
         --bs-btn-active-bg: #5a23c8;
         --bs-btn-active-border-color: #5a23c8;
         }
         .bd-mode-toggle {
         z-index: 1500;
         }
         .bd-mode-toggle .bi {
         width: 1em;
         height: 1em;
         }
         .bd-mode-toggle .dropdown-menu .active .bi {
         display: block !important;
         }
		 .cardtile-header{
			font-size: 22px;
		 }
		 .cardtile-header2{
			font-size: 16px;
		 }
      </style>
   </head>
   <body>
      <svg xmlns="http://www.w3.org/2000/svg" class="d-none">
         <symbol id="check2" viewBox="0 0 16 16">
            <path d="M13.854 3.646a.5.5 0 0 1 0 .708l-7 7a.5.5 0 0 1-.708 0l-3.5-3.5a.5.5 0 1 1 .708-.708L6.5 10.293l6.646-6.647a.5.5 0 0 1 .708 0z"></path>
         </symbol>
         <symbol id="circle-half" viewBox="0 0 16 16">
            <path d="M8 15A7 7 0 1 0 8 1v14zm0 1A8 8 0 1 1 8 0a8 8 0 0 1 0 16z"></path>
         </symbol>
         <symbol id="moon-stars-fill" viewBox="0 0 16 16">
            <path
               d="M6 .278a.768.768 0 0 1 .08.858 7.208 7.208 0 0 0-.878 3.46c0 4.021 3.278 7.277 7.318 7.277.527 0 1.04-.055 1.533-.16a.787.787 0 0 1 .81.316.733.733 0 0 1-.031.893A8.349 8.349 0 0 1 8.344 16C3.734 16 0 12.286 0 7.71 0 4.266 2.114 1.312 5.124.06A.752.752 0 0 1 6 .278z"
               ></path>
            <path
               d="M10.794 3.148a.217.217 0 0 1 .412 0l.387 1.162c.173.518.579.924 1.097 1.097l1.162.387a.217.217 0 0 1 0 .412l-1.162.387a1.734 1.734 0 0 0-1.097 1.097l-.387 1.162a.217.217 0 0 1-.412 0l-.387-1.162A1.734 1.734 0 0 0 9.31 6.593l-1.162-.387a.217.217 0 0 1 0-.412l1.162-.387a1.734 1.734 0 0 0 1.097-1.097l.387-1.162zM13.863.099a.145.145 0 0 1 .274 0l.258.774c.115.346.386.617.732.732l.774.258a.145.145 0 0 1 0 .274l-.774.258a1.156 1.156 0 0 0-.732.732l-.258.774a.145.145 0 0 1-.274 0l-.258-.774a1.156 1.156 0 0 0-.732-.732l-.774-.258a.145.145 0 0 1 0-.274l.774-.258c.346-.115.617-.386.732-.732L13.863.1z"
               ></path>
         </symbol>
         <symbol id="sun-fill" viewBox="0 0 16 16">
            <path
               d="M8 12a4 4 0 1 0 0-8 4 4 0 0 0 0 8zM8 0a.5.5 0 0 1 .5.5v2a.5.5 0 0 1-1 0v-2A.5.5 0 0 1 8 0zm0 13a.5.5 0 0 1 .5.5v2a.5.5 0 0 1-1 0v-2A.5.5 0 0 1 8 13zm8-5a.5.5 0 0 1-.5.5h-2a.5.5 0 0 1 0-1h2a.5.5 0 0 1 .5.5zM3 8a.5.5 0 0 1-.5.5h-2a.5.5 0 0 1 0-1h2A.5.5 0 0 1 3 8zm10.657-5.657a.5.5 0 0 1 0 .707l-1.414 1.415a.5.5 0 1 1-.707-.708l1.414-1.414a.5.5 0 0 1 .707 0zm-9.193 9.193a.5.5 0 0 1 0 .707L3.05 13.657a.5.5 0 0 1-.707-.707l1.414-1.414a.5.5 0 0 1 .707 0zm9.193 2.121a.5.5 0 0 1-.707 0l-1.414-1.414a.5.5 0 0 1 .707-.707l1.414 1.414a.5.5 0 0 1 0 .707zM4.464 4.465a.5.5 0 0 1-.707 0L2.343 3.05a.5.5 0 1 1 .707-.707l1.414 1.414a.5.5 0 0 1 0 .708z"
               ></path>
         </symbol>
      </svg>
      <div
         class="dropdown position-fixed bottom-0 end-0 mb-3 me-3 bd-mode-toggle"
         >
         <button
            class="btn btn-bd-primary py-2 dropdown-toggle d-flex align-items-center"
            id="bd-theme"
            type="button"
            aria-expanded="false"
            data-bs-toggle="dropdown"
            aria-label="Toggle theme (auto)"
            >
            <svg class="bi my-1 theme-icon-active" aria-hidden="true">
               <use href="#circle-half"></use>
            </svg>
            <span class="visually-hidden" id="bd-theme-text">Toggle theme</span>
         </button>
         <ul
            class="dropdown-menu dropdown-menu-end shadow"
            aria-labelledby="bd-theme-text"
            >
            <li>
               <button
                  type="button"
                  class="dropdown-item d-flex align-items-center"
                  data-bs-theme-value="light"
                  aria-pressed="false"
                  >
                  <svg class="bi me-2 opacity-50" aria-hidden="true">
                     <use href="#sun-fill"></use>
                  </svg>
                  Light
                  <svg class="bi ms-auto d-none" aria-hidden="true">
                     <use href="#check2"></use>
                  </svg>
               </button>
            </li>
            <li>
               <button
                  type="button"
                  class="dropdown-item d-flex align-items-center"
                  data-bs-theme-value="dark"
                  aria-pressed="false"
                  >
                  <svg class="bi me-2 opacity-50" aria-hidden="true">
                     <use href="#moon-stars-fill"></use>
                  </svg>
                  Dark
                  <svg class="bi ms-auto d-none" aria-hidden="true">
                     <use href="#check2"></use>
                  </svg>
               </button>
            </li>
            <li>
               <button
                  type="button"
                  class="dropdown-item d-flex align-items-center active"
                  data-bs-theme-value="auto"
                  aria-pressed="true"
                  >
                  <svg class="bi me-2 opacity-50" aria-hidden="true">
                     <use href="#circle-half"></use>
                  </svg>
                  Auto
                  <svg class="bi ms-auto d-none" aria-hidden="true">
                     <use href="#check2"></use>
                  </svg>
               </button>
            </li>
         </ul>
      </div>
      <header data-bs-theme="dark">
         <div class="collapse text-bg-dark" id="navbarHeader">
            <div class="container">
               <div class="row">
                  <div class="col-sm-8 col-md-7 py-4">
                     <h4>About</h4>
                     <p class="text-body-secondary">
						<br />
						<b>Tropism Analysis Package: Interactive Machine Learning Software to Identify Viral Host Factors Through Single-Cell Host-Virus mRNA Profiling</b><br />
						Authors: <i>Kenny Pavan, Lamya Ben Ameur, Zach Goode, Emily Tiedemann, Elizabeth Kolb, Arpiar Saunders</i><br />
                     </p>
                  </div>
                  <div class="col-sm-4 offset-md-1 py-4">
                     <h4>Contact</h4>
                     <ul class="list-unstyled">
                        <li><a href="https://github.com/ArpiarSaundersLab/tap" target="_blank" class="text-white">Github</a></li>
                     </ul>
                  </div>
               </div>
            </div>
         </div>
         <div class="navbar navbar-dark bg-dark shadow-sm">
            <div class="container">
				<a href="#" class="navbar-brand d-flex align-items-center">
					<svg
						xmlns="http://www.w3.org/2000/svg"
						width="20"
						height="20"
						fill="none"
						stroke="currentColor"
						stroke-linecap="round"
						stroke-linejoin="round"
						stroke-width="2"
						viewBox="0 0 24 24"
						class="me-2">
						<polygon points="12,2 20,7 20,17 12,22 4,17 4,7" />
						<path d="M12 2 L12 22" />
						<path d="M4 7 L20 17" />
						<path d="M20 7 L4 17" />
						<circle cx="12" cy="12" r="2" />
					</svg>
					<strong>AAVdb</strong>
			    </a>
               <button
					class="navbar-toggler"
					type="button"
					data-bs-toggle="collapse"
					data-bs-target="#navbarHeader"
					aria-controls="navbarHeader"
					aria-expanded="false" aria-label="Toggle navigation">
               		<span class="navbar-toggler-icon"></span>
               </button>
            </div>
         </div>
      </header>
      <main>
         <section class="py-5 text-center container">
            <div class="row py-lg-5">
               <div class="col-lg-7 col-md-8 mx-auto">
                  <h1 class="fw-light">Interactive Single-Cell  Viral Tropism Portal</h1>
                  <p class="lead text-body-secondary">
                     Below are some example AAV tropism datasets generated using the Tropism Analysis Package (TAP). 
					 Please explore the datasets and feel free to reach out with any questions or suggestions!
                  </p>
				  <a href="tap_manuscript.pdf" class="btn btn-lg btn-primary my-2" target="_blank">Manuscript</a>
				  <a href="tap_methods.pdf" class="btn btn-lg btn-secondary my-2 ms-3" target="_blank">Methods</a>
               </div>
            </div>
         </section>
         <div class="album py-5 bg-body-tertiary">
            <div class="container">
               <div class="row row-cols-1 row-cols-sm-2 row-cols-md-3 g-3">
                  <div class="col">
                     <div class="card shadow-sm">
						<a href="taps/1b.html" target="_blank">
                        <svg
                           aria-label="Placeholder: Thumbnail"
                           class="bd-placeholder-img card-img-top"
                           height="225"
                           preserveAspectRatio="xMidYMid slice"
                           role="img"
                           width="100%"
                           xmlns="http://www.w3.org/2000/svg">
                           <title>In-Vitro AAVs</title>
                           <rect width="100%" height="100%" fill="#16526F"></rect>
                           <text x="50%" y="50%" fill="#eceeef" dy=".3em" class="cardtile-header">In-Vitro AAVs</text>
                        </svg>
						</a>
                        <div class="card-body">
                           <p class="card-text">
                              A collection of 11 naturally occurring and engineered AAVs were used to transduce primary mixed mouse cortical cultures.
                           </p>
							<div class="btn-group">
								<a href="taps/1b.html" target="_blank" class="btn btn-sm btn-outline-secondary">Single</a>
								<a href="taps/1a.html" target="_blank" class="btn btn-sm btn-outline-secondary ms-2">Replicates</a>
							</div>
                           <div class="text-end">
                              <small class="text-body-secondary">6,968 | 16,418 cells</small>
                           </div>
                        </div>
                     </div>
                  </div>
                  <div class="col">
                     <div class="card shadow-sm">
						<a href="taps/2.html" target="_blank">
                        <svg
                           aria-label="Placeholder: Thumbnail"
                           class="bd-placeholder-img card-img-top"
                           height="225"
                           preserveAspectRatio="xMidYMid slice"
                           role="img"
                           width="100%"
                           xmlns="http://www.w3.org/2000/svg">
                           <title>In-Vivo AAVs</title>
                           <rect width="100%" height="100%" fill="#178D17"></rect>
                           <text x="50%" y="50%" fill="#eceeef" dy=".3em" class="cardtile-header">In-Vivo AAVs</text>
                        </svg>
						</a>
                        <div class="card-body">
                           <p class="card-text">
                              14 barcoded AAVs introduced into the lateral ventricle of E13.5 mouse embryos. Samples were pooled and single-cell RNA-seq was performed.
                           </p>
                           <div class="text-end">
                              <small class="text-body-secondary">12,875 cells</small>
                           </div>
                        </div>
                     </div>
                  </div>
                  <div class="col">
                     <div class="card shadow-sm">
						<a href="taps/3.html" target="_blank">
                        <svg
                           aria-label="Placeholder: Thumbnail"
                           class="bd-placeholder-img card-img-top"
                           height="225"
                           preserveAspectRatio="xMidYMid slice"
                           role="img"
                           width="100%"
                           xmlns="http://www.w3.org/2000/svg">
                           <title>Extended Use Case</title>
                           <rect width="100%" height="100%" fill="#B06C1C"></rect>
                           <text x="50%" y="50%" fill="#eceeef" dy=".3em" class="cardtile-header">Extended Use Case</text>
						   <text x="50%" y="65%" fill="#eceeef" dy=".2em" class="cardtile-header2">SAD-B19ΔG</text>
                        </svg>
						</a>
                        <div class="card-body">
                           <p class="card-text">
                              Rabies Virus SAD-B19ΔG was used to transduce primary mouse cortical cultures, then single-cell RNA-seq was performed and run through TAP.
                           </p>
							<div class="btn-group">
								<a href="taps/3.html" target="_blank" class="btn btn-sm btn-outline-secondary">Single</a>
								<a href="taps/3a.html" target="_blank" class="btn btn-sm btn-outline-secondary ms-2">Replicates</a>
							</div>
                           <div class="text-end">
                              <small class="text-body-secondary">12,680 cells</small>
                           </div>
                        </div>
                     </div>
                  </div>

               </div>
            </div>
         </div>
      </main>
      <footer class="text-body-secondary py-5">
         <div class="container">
            <p class="mb-1">AAVdb.com</p>
            <p class="mb-0">
              Created by the <a href="https://www.arpiarsaunderslab.org/" target="_blank">Saunders Lab</a> in the <a href="https://www.ohsu.edu/vollum-institute" target="_blank">Vollum Institute</a> @ <a href="https://www.ohsu.edu/" target="_blank">OHSU</a> <br />
			  
            </p>
         </div>
      </footer>
      <script
         src="assets/dist/js/bootstrap.bundle.min.js"
         class="astro-vvvwv3sm"
         ></script>
   </body>
</html>