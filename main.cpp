#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h" // https://github.com/nothings/stb/blob/master/stb_image_write.h

// Explicitly include all necessary SFML headers for clarity
// SFML 2.5.1 uses sf::Uint8 directly.
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>   // For sf::Event, sf::VideoMode
#include <SFML/System.hpp>   // For sf::Vector2, sf::Clock, sf::Time, sf::Uint8

#include <imgui.h>
#include <imgui-SFML.h>

#include <vector>
#include <cmath>
#include <iostream>
#include <thread>
#include <mutex>
#include <atomic>
#include <iomanip>     // For std::setprecision
#include <cstring>     // For std::memcpy
#include <algorithm>   // For std::min and std::max

// Define M_PI if not available (e.g., in some C++ versions or environments)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Constants (Geometric Units: G=c=1) ---
const double M_DEFAULT = 1.0; // Default Black Hole Mass
const double R_SCHWARZSCHILD_FACTOR = 2.0; // r_s = 2M
const double R_ISCO_FACTOR = 6.0;         // r_isco = 6M
const double R_OUTER_DISK_FACTOR = 100.0; // Outer radius of the accretion disk (e.g., 100M)
const double EPSILON = 1e-9; // Small value to prevent division by zero or log(0)
const double RK4_STEP_SIZE_DEFAULT = 0.05; // Default integration step size for RK4
const int MAX_RK4_STEPS_DEFAULT = 5000;    // Default max steps per ray to prevent infinite loops

// --- Global UI State ---
int render_image_width = 720;
int render_image_height = 720;
int gui_panel_width = 300; // Width of the GUI control panel

float current_camera_angle_deg = 90.0f; // Initial Camera Angle
float current_bh_mass = M_DEFAULT;      // Initial Black Hole Mass

std::atomic<bool> render_in_progress = false;
std::atomic<bool> render_cancel_flag = false; // Flag to signal rendering thread to stop

// Use sf::Uint8 for pixel data in SFML 2.5.1
std::vector<sf::Uint8> pixel_buffer_rgba; // Stores RGBA pixel data
sf::Texture rendered_texture; // Global texture object
sf::Sprite rendered_sprite;  // Global sprite object
std::mutex pixel_buffer_mutex; // Protects pixel_buffer_rgba and rendered_texture updates

// --- Physics Data Structures ---

// Represents a photon's state in 8 dimensions (position and momentum)
struct PhotonState {
    double t, r, theta, phi;         // Position in Schwarzschild coordinates
    double pt, pr, ptheta, pphi;     // Covariant momenta conjugate to position
};

// --- Schwarzschild Metric Functions ---

// Get the g_tt component of the Schwarzschild metric
double g_tt(double r, double M) {
    if (r <= R_SCHWARZSCHILD_FACTOR * M + EPSILON) return -EPSILON; // Avoid singularity
    return -(1.0 - R_SCHWARZSCHILD_FACTOR * M / r);
}

// Get the g_rr component
double g_rr(double r, double M) {
    if (r <= R_SCHWARZSCHILD_FACTOR * M + EPSILON) return 1.0 / EPSILON; // Avoid singularity
    return 1.0 / (1.0 - R_SCHWARZSCHILD_FACTOR * M / r);
}

// Get the g_theta_theta component
double g_theta_theta(double r) {
    return r * r;
}

// Get the g_phi_phi component
double g_phi_phi(double r, double theta) {
    double sin_theta = std::sin(theta);
    if (std::abs(sin_theta) < EPSILON) {
        return r * r * EPSILON; // Prevent zero or near-zero results at poles
    }
    return r * r * sin_theta * sin_theta;
}

// Get the inverse metric component G^tt
double G_tt(double r, double M) {
    double val_g_tt = g_tt(r, M);
    if (std::abs(val_g_tt) < EPSILON) return 1.0 / -EPSILON; // Ensure non-zero inverse
    return 1.0 / val_g_tt;
}

// Get the inverse metric component G^rr
double G_rr(double r, double M) {
    double val_g_rr = g_rr(r, M);
    if (std::abs(val_g_rr) < EPSILON) return 1.0 / EPSILON; // Ensure non-zero inverse
    return 1.0 / val_g_rr;
}

// Get the inverse metric component G^theta_theta
double G_theta_theta(double r) {
    double val_g_theta_theta = g_theta_theta(r);
    if (std::abs(val_g_theta_theta) < EPSILON) return 1.0 / EPSILON; // Ensure non-zero inverse
    return 1.0 / val_g_theta_theta;
}

// Get the inverse metric component G^phi_phi
double G_phi_phi(double r, double theta) {
    double val_g_phi_phi = g_phi_phi(r, theta);
    if (std::abs(val_g_phi_phi) < EPSILON) return 1.0 / EPSILON; // Ensure non-zero inverse
    return 1.0 / val_g_phi_phi;
}

// --- Christoffel Symbols (non-zero ones for Schwarzschild) ---
double Gamma_r_tt(double r, double M) {
    if (r < EPSILON) return 0.0;
    return (M / (r * r)) * (1.0 - R_SCHWARZSCHILD_FACTOR * M / r);
}

double Gamma_t_tr(double r, double M) {
    if (r < EPSILON || std::abs(1.0 - R_SCHWARZSCHILD_FACTOR * M / r) < EPSILON) return 0.0;
    return M / (r * r * (1.0 - R_SCHWARZSCHILD_FACTOR * M / r));
}

double Gamma_r_rr(double r, double M) {
    if (r < EPSILON || std::abs(1.0 - R_SCHWARZSCHILD_FACTOR * M / r) < EPSILON) return 0.0;
    return -M / (r * r * (1.0 - R_SCHWARZSCHILD_FACTOR * M / r));
}

double Gamma_r_theta_theta(double r) {
    return -r;
}

double Gamma_r_phi_phi(double r, double theta) {
    return -r * std::sin(theta) * std::sin(theta);
}

double Gamma_theta_r_theta(double r) {
    if (r < EPSILON) return 0.0;
    return 1.0 / r;
}

double Gamma_theta_phi_phi(double theta) {
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);
    if (std::abs(sin_theta) < EPSILON) return 0.0; // Avoid division by zero from tan() later if cos != 0
    return -sin_theta * cos_theta;
}

double Gamma_phi_r_phi(double r) {
    if (r < EPSILON) return 0.0;
    return 1.0 / r;
}

double Gamma_phi_theta_phi(double theta) {
    double tan_theta = std::tan(theta);
    if (std::abs(tan_theta) < EPSILON) return 0.0; // Avoid division by zero
    return 1.0 / tan_theta;
}

// --- Geodesic Equations ---
PhotonState computeDerivatives(const PhotonState& s, double M) {
    PhotonState ds;

    // Clamp theta to prevent issues at poles (theta=0 or theta=PI) for metric components
    double safe_theta = std::max(EPSILON, std::min(s.theta, M_PI - EPSILON));

    // Calculate contravariant momentum components p^mu from covariant p_mu
    double pt_contra = G_tt(s.r, M) * s.pt;
    double pr_contra = G_rr(s.r, M) * s.pr;
    double ptheta_contra = G_theta_theta(s.r) * s.ptheta;
    double pphi_contra = G_phi_phi(s.r, safe_theta) * s.pphi;

    // dx^mu / d_lambda = p^mu
    ds.t = pt_contra;
    ds.r = pr_contra;
    ds.theta = ptheta_contra;
    ds.phi = pphi_contra;

    // dp_mu / d_lambda = -Gamma^alpha_mu_beta * p^alpha * p^beta
    ds.pt = -(
        Gamma_t_tr(s.r, M) * pt_contra * pr_contra * 2.0 // 2 * Gamma^t_tr * p^t * p^r
        );

    ds.pr = -(
        Gamma_r_tt(s.r, M) * pt_contra * pt_contra +
        Gamma_r_rr(s.r, M) * pr_contra * pr_contra +
        Gamma_r_theta_theta(s.r) * ptheta_contra * ptheta_contra +
        Gamma_r_phi_phi(s.r, safe_theta) * pphi_contra * pphi_contra
        );

    ds.ptheta = -(
        Gamma_theta_r_theta(s.r) * ptheta_contra * pr_contra * 2.0 +
        Gamma_theta_phi_phi(safe_theta) * pphi_contra * pphi_contra
        );

    ds.pphi = -(
        Gamma_phi_r_phi(s.r) * pphi_contra * pr_contra * 2.0 +
        Gamma_phi_theta_phi(safe_theta) * pphi_contra * ptheta_contra * 2.0
        );

    return ds;
}

PhotonState operator+(const PhotonState& a, const PhotonState& b) {
    return {
        a.t + b.t, a.r + b.r, a.theta + b.theta, a.phi + b.phi,
        a.pt + b.pt, a.pr + b.pr, a.ptheta + b.ptheta, a.pphi + b.pphi
    };
}

PhotonState operator*(const PhotonState& s, double scalar) {
    return {
        s.t * scalar, s.r * scalar, s.theta * scalar, s.phi * scalar,
        s.pt * scalar, s.pr * scalar, s.ptheta * scalar, s.pphi * scalar
    };
}

void rungeKutta4(PhotonState& s, double h, double M) {
    PhotonState k1 = computeDerivatives(s, M);
    PhotonState k2 = computeDerivatives(s + k1 * (h / 2.0), M);
    PhotonState k3 = computeDerivatives(s + k2 * (h / 2.0), M);
    PhotonState k4 = computeDerivatives(s + k3 * h, M);

    s = s + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);

    // Normalize phi to [0, 2PI)
    s.phi = std::fmod(s.phi, 2.0 * M_PI);
    if (s.phi < 0) s.phi += 2.0 * M_PI;

    // Normalize theta to [0, PI] (reflecting across equatorial plane if needed)
    s.theta = std::fmod(s.theta, 2.0 * M_PI); // First bring into [0, 2PI)
    if (s.theta < 0) s.theta += 2.0 * M_PI;
    if (s.theta > M_PI) { // If it's in (PI, 2PI), reflect it
        s.theta = 2.0 * M_PI - s.theta;
    }
    // Clamp theta to avoid numerical issues very close to poles
    s.theta = std::max(EPSILON, std::min(s.theta, M_PI - EPSILON));
}

// --- Ray Tracing Core ---
struct RayHitInfo {
    bool is_captured;
    bool is_disk_hit;
    double r_disk;
};

RayHitInfo traceRay(double camera_theta_0_rad, double pixel_x_norm, double pixel_y_norm,
    double M, double r_outer_disk, double rk4_step_size, int max_rk4_steps,
    std::atomic<bool>& cancel_flag) {
    RayHitInfo hit = { false, false, 0.0 };

    double obs_r = r_outer_disk * 20.0; // Observer far away
    double obs_theta = camera_theta_0_rad; // Observer's inclination
    double obs_phi = M_PI; // Arbitrary observer phi (looking "backwards")

    // Ensure observer theta is not exactly 0 or PI to avoid division by sin(theta) later
    obs_theta = std::max(EPSILON, std::min(obs_theta, M_PI - EPSILON));

    PhotonState s = {
        0.0,        // t (arbitrary, start at 0)
        obs_r,      // r
        obs_theta,  // theta
        obs_phi,    // phi
        0.0, 0.0, 0.0, 0.0 // Placeholder momenta
    };

    s.pt = -1.0; // Assume unit energy for the photon, moving backwards in time (p_t = -E)

    // Initial momenta based on screen coordinates (impact parameters)
    // These values represent the initial covariant angular momenta.
    s.pphi = pixel_x_norm * obs_r * std::sin(obs_theta);
    s.ptheta = pixel_y_norm * obs_r;

    // Calculate initial p_r from the null geodesic condition (g^mu nu p_mu p_nu = 0)
    // g_tt (p_t)^2 + g_rr (p_r)^2 + g_theta_theta (p_theta)^2 + g_phi_phi (p_phi)^2 = 0
    // Solving for (p_r)^2:
    // (p_r)^2 = -(g_tt (p_t)^2 + g_theta_theta (p_theta)^2 + g_phi_phi (p_phi)^2) / g_rr

    double val_g_tt_at_s_r = g_tt(s.r, M);
    double val_g_rr_at_s_r = g_rr(s.r, M);
    double val_g_theta_theta_at_s_r = g_theta_theta(s.r);
    double val_g_phi_phi_at_s_r_theta = g_phi_phi(s.r, s.theta);

    double p_r_squared_numerator = -(val_g_tt_at_s_r * s.pt * s.pt +
        val_g_theta_theta_at_s_r * s.ptheta * s.ptheta +
        val_g_phi_phi_at_s_r_theta * s.pphi * s.pphi);

    if (val_g_rr_at_s_r < EPSILON || p_r_squared_numerator < 0) {
        // This ray is problematic (e.g., trying to start inside event horizon, or invalid momenta)
        // Treat as captured or escaped, or give it a minimal inward push.
        s.pr = -std::sqrt(EPSILON); // Small inward momentum to attempt continuation
    }
    else {
        s.pr = -std::sqrt(p_r_squared_numerator / val_g_rr_at_s_r); // Negative for an inward-going ray
    }

    for (int step = 0; step < max_rk4_steps; ++step) {
        if (cancel_flag.load()) return hit;

        // Check for capture by black hole (r <= Schwarzschild radius)
        if (s.r <= R_SCHWARZSCHILD_FACTOR * M) {
            hit.is_captured = true;
            return hit;
        }
        // Also if r somehow becomes negative or effectively zero due to numerical instability
        if (s.r < EPSILON) {
            hit.is_captured = true;
            return hit;
        }

        // Check for hitting the accretion disk (r >= ISCO and near equatorial plane)
        double angular_tolerance = 5.0 * M_PI / 180.0; // 5 degrees in radians for hitting disk
        if (s.r >= R_ISCO_FACTOR * M && s.r <= r_outer_disk &&
            std::abs(s.theta - M_PI / 2.0) < angular_tolerance) {
            hit.is_disk_hit = true;
            hit.r_disk = s.r;
            return hit;
        }

        // If ray goes too far out beyond observer's assumed max distance, assume it missed
        if (s.r > obs_r * 1.5) {
            return hit;
        }

        rungeKutta4(s, rk4_step_size, M);
    }

    // If max steps reached without hitting anything (ray went to infinity or got stuck)
    return hit;
}

// --- Luminosity Function (Simplified) ---
double getDiskLuminosity(double r_disk, double M) {
    if (r_disk < R_ISCO_FACTOR * M) { // Below ISCO, no emission
        return 0.0;
    }

    double r_norm = r_disk / M; // Normalize r by M as in the formula

    // Simplified Page & Thorne formula (Luminet Eq. 16, without M_dot constant and log term)
    double fs_value = (r_norm - 3.0) / (r_norm * r_norm * std::sqrt(r_norm) + EPSILON);

    return std::max(0.0, fs_value); // Luminosity cannot be negative
}

// Helper functions for float to byte array conversion (due to storing luminosity in pixel buffer)
void float_to_byte_array(float f, sf::Uint8* bytes) {
    static_assert(sizeof(float) == 4, "float is not 4 bytes!"); // Ensure float is 4 bytes
    std::memcpy(bytes, &f, sizeof(float));
}

float byte_array_to_float(const sf::Uint8* bytes) {
    float f;
    static_assert(sizeof(float) == 4, "float is not 4 bytes!"); // Ensure float is 4 bytes
    std::memcpy(&f, bytes, sizeof(float));
    return f;
}

// --- Rendering Thread Function ---
void render_thread_function(float angle_deg, float bh_mass, int width, int height) {
    // Local buffer for this thread's rendering data
    std::vector<sf::Uint8> local_pixel_buffer(static_cast<size_t>(width) * height * 4);

    double camera_angle_rad = angle_deg * M_PI / 180.0; // Manual radians conversion for C++17/SFML 2.5.1
    double r_outer_disk_val = R_OUTER_DISK_FACTOR * bh_mass;
    double max_luminosity = 0.0; // To find max flux for normalization

    // First pass: trace rays and find max luminosity for normalization
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (render_cancel_flag.load()) { // Check cancel flag
                std::fill(local_pixel_buffer.begin(), local_pixel_buffer.end(), 0);
                render_in_progress.store(false);
                render_cancel_flag.store(false);
                return;
            }

            // Normalize pixel coordinates to [-1, 1] range.
            // Using (width - 1.0) and (height - 1.0) for accurate mapping from 0 to max pixel.
            double normalized_x = (static_cast<double>(x) / (width - 1.0)) * 2.0 - 1.0;
            double normalized_y = (static_cast<double>(y) / (height - 1.0)) * 2.0 - 1.0;

            RayHitInfo hit = traceRay(camera_angle_rad, normalized_x, normalized_y,
                bh_mass, r_outer_disk_val, RK4_STEP_SIZE_DEFAULT, MAX_RK4_STEPS_DEFAULT,
                render_cancel_flag);

            double luminosity = 0.0;
            if (hit.is_disk_hit) {
                luminosity = getDiskLuminosity(hit.r_disk, bh_mass);
                luminosity = std::max(0.0, luminosity); // Ensure non-negative
                if (luminosity > max_luminosity) {
                    max_luminosity = luminosity;
                }
            }

            // Store raw luminosity as float in the RGBA buffer for the second pass
            // We use the first 4 bytes of each pixel to store a float.
            size_t index_byte = (static_cast<size_t>(y) * width + x) * 4;
            if (index_byte + 3 < local_pixel_buffer.size()) { // Bounds check
                float_to_byte_array(static_cast<float>(luminosity), &local_pixel_buffer[index_byte]);
            }
        }
    }

    // Second pass: apply normalization and convert to grayscale
    if (max_luminosity < EPSILON) max_luminosity = 1.0; // Avoid division by zero

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (render_cancel_flag.load()) { // Check cancel flag again
                std::fill(local_pixel_buffer.begin(), local_pixel_buffer.end(), 0); // Clear on cancel
                render_in_progress.store(false);
                render_cancel_flag.store(false);
                return;
            }

            size_t index_byte = (static_cast<size_t>(y) * width + x) * 4;
            if (index_byte + 3 < local_pixel_buffer.size()) { // Bounds check
                float stored_luminosity = byte_array_to_float(&local_pixel_buffer[index_byte]);

                sf::Uint8 gray_value;
                if (stored_luminosity == 0.0) { // Black hole, empty space, or no hit
                    gray_value = 0;
                }
                else {
                    double normalized_brightness = stored_luminosity / max_luminosity;
                    // Apply a simple gamma correction for better visual appearance (e.g., gamma = 2.2)
                    normalized_brightness = std::pow(normalized_brightness, 1.0 / 2.2);
                    gray_value = static_cast<sf::Uint8>(std::min(normalized_brightness * 255.0, 255.0));
                }

                local_pixel_buffer[index_byte + 0] = gray_value;
                local_pixel_buffer[index_byte + 1] = gray_value;
                local_pixel_buffer[index_byte + 2] = gray_value;
                local_pixel_buffer[index_byte + 3] = 255; // Alpha
            }
        }
    }

    // --- Update shared pixel buffer and texture ---
    std::lock_guard<std::mutex> lock(pixel_buffer_mutex);
    pixel_buffer_rgba = std::move(local_pixel_buffer); // Move data efficiently

    // SFML 2.5.1: sf::Texture::create takes unsigned int width, height
    if (rendered_texture.getSize().x != static_cast<unsigned int>(width) || rendered_texture.getSize().y != static_cast<unsigned int>(height)) {
        std::cerr << "Texture size mismatch or not created, recreating texture." << std::endl;
        rendered_texture.create(static_cast<unsigned int>(width), static_cast<unsigned int>(height));
        rendered_sprite.setTexture(rendered_texture, true); // Re-attach texture to sprite
    }
    rendered_texture.update(pixel_buffer_rgba.data());

    render_in_progress.store(false); // Rendering finished
    render_cancel_flag.store(false); // Reset cancel flag
    std::cout << "Rendering finished." << std::endl;
}

// --- Image Saving Function ---
void save_image_to_file(const std::string& filename, int width, int height, const std::vector<sf::Uint8>& buffer) {
    if (buffer.empty()) {
        std::cerr << "No image to save." << std::endl;
        return;
    }

    // stbi_write_png takes (filename, width, height, num_channels, data, stride_in_bytes)
    // For RGBA (4 channels), stride is width * 4.
    int result = stbi_write_png(filename.c_str(), width, height, 4, buffer.data(), width * 4);

    if (result) {
        std::cout << "Image saved successfully to " << filename << std::endl;
    }
    else {
        std::cerr << "Failed to save image to " << filename << std::endl;
    }
}

// --- Main Application Loop ---
int main() {
    // SFML 2.6.1: sf::VideoMode constructor takes width, height directly.
    sf::RenderWindow window(sf::VideoMode(static_cast<unsigned int>(render_image_width + gui_panel_width), static_cast<unsigned int>(render_image_height)), "BLACK HOLE COMPUTATION");
    window.setFramerateLimit(60); // Limit to 60 FPS

    // Initialize ImGui-SFML
    if (!ImGui::SFML::Init(window)) {
        std::cerr << "Failed to initialize ImGui-SFML! Check library versions." << std::endl;
        return 1;
    }

    // Initialize pixel buffer and texture immediately within main.
    // SFML 2.6.1: sf::Texture::create takes unsigned int width, height.
    pixel_buffer_rgba.resize(static_cast<size_t>(render_image_width) * render_image_height * 4);
    rendered_texture.create(static_cast<unsigned int>(render_image_width), static_cast<unsigned int>(render_image_height));

    // sf::Sprite needs a texture, set it immediately after texture is created.
    rendered_sprite.setTexture(rendered_texture, true); // true to reset texture rect

    // Start initial render
    render_in_progress.store(true); // Mark as rendering in progress
    std::thread initial_render_worker(render_thread_function, current_camera_angle_deg, current_bh_mass, render_image_width, render_image_height);
    initial_render_worker.detach(); // Let the thread run independently

    sf::Clock deltaClock; // For ImGui-SFML update

    while (window.isOpen()) {
        sf::Event event; // SFML 2.6.1 Event is default constructible and pollEvent takes sf::Event&
        while (window.pollEvent(event)) {
            ImGui::SFML::ProcessEvent(event); // ImGui-SFML should also take sf::Event&
            if (event.type == sf::Event::Closed)
                window.close();
        }

        ImGui::SFML::Update(window, deltaClock.restart());

        // --- GUI Layout ---
        ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f)); // Place GUI panel at top-left
        // Ensure ImVec2 takes floats
        ImGui::SetNextWindowSize(ImVec2(static_cast<float>(gui_panel_width), static_cast<float>(render_image_height))); // Fixed size
        ImGui::Begin("Controls", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);

        bool rendering = render_in_progress.load(); // Get current rendering state

        // Top-right Copy/Save button
        float button_width = 80.0f;
        float button_height = 25.0f;
        ImGui::SetCursorPosX(ImGui::GetWindowWidth() - button_width - ImGui::GetStyle().WindowPadding.x);
        ImGui::SetCursorPosY(ImGui::GetStyle().WindowPadding.y);
        // Disable if rendering, or if buffer is empty/all black (initial state)
        bool buffer_is_empty_or_black = pixel_buffer_rgba.empty() || (pixel_buffer_rgba.size() >= 4 && pixel_buffer_rgba[0] == 0 && pixel_buffer_rgba[1] == 0 && pixel_buffer_rgba[2] == 0);
        ImGui::BeginDisabled(rendering || buffer_is_empty_or_black);
        if (ImGui::Button("Copy/Save", ImVec2(button_width, button_height))) {
            std::lock_guard<std::mutex> lock(pixel_buffer_mutex); // Lock to save buffer safely
            save_image_to_file("rendered_black_hole.png", render_image_width, render_image_height, pixel_buffer_rgba);
        }
        ImGui::EndDisabled();

        ImGui::Separator(); // Visual separator

        // Render and Cancel Buttons
        ImGui::SetCursorPosY(ImGui::GetCursorPosY() + 10.0f); // Some spacing
        ImGui::BeginDisabled(rendering); // Disable Render button if rendering
        // Ensure float arguments for ImVec2
        if (ImGui::Button("Render", ImVec2(ImGui::GetContentRegionAvail().x / 2.0f - 5.0f, 0.0f))) {
            // Start rendering in a new thread
            render_cancel_flag.store(false); // Reset cancel flag
            render_in_progress.store(true);  // Mark as rendering in progress
            std::thread render_worker(render_thread_function, current_camera_angle_deg, current_bh_mass, render_image_width, render_image_height);
            render_worker.detach();
            std::cout << "Rendering started with angle: " << current_camera_angle_deg << ", mass: " << current_bh_mass << std::endl;
        }
        ImGui::EndDisabled();
        ImGui::SameLine();
        ImGui::BeginDisabled(!rendering); // Disable Cancel button if not rendering
        // Ensure float arguments for ImVec2
        if (ImGui::Button("Cancel", ImVec2(ImGui::GetContentRegionAvail().x / 2.0f - 5.0f, 0.0f))) {
            render_cancel_flag.store(true); // Signal render thread to stop
            std::cout << "Rendering cancelled." << std::endl;
        }
        ImGui::EndDisabled();

        ImGui::Separator(); // Visual separator

        // Theta Slider and Input
        ImGui::Text("Theta:");
        ImGui::BeginDisabled(rendering); // Disable if rendering
        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x - 60.0f); // Make slider wide, leave space for text input
        ImGui::SliderFloat("##ThetaSlider", &current_camera_angle_deg, 0.0f, 180.0f, "%.1f deg");
        ImGui::SameLine();
        ImGui::PopItemWidth();
        ImGui::PushItemWidth(50.0f); // Small width for text input
        ImGui::InputFloat("##ThetaInput", &current_camera_angle_deg, 0.0f, 0.0f, "%.1f");
        ImGui::PopItemWidth();
        ImGui::EndDisabled();

        // Mass Slider and Input
        ImGui::Text("Mass:");
        ImGui::BeginDisabled(rendering); // Disable if rendering
        ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x - 60.0f);
        ImGui::SliderFloat("##MassSlider", &current_bh_mass, 0.5f, 5.0f, "%.2f M"); // Adjust min/max as needed
        ImGui::SameLine();
        ImGui::PopItemWidth();
        ImGui::PushItemWidth(50.0f);
        ImGui::InputFloat("##MassInput", &current_bh_mass, 0.0f, 0.0f, "%.2f");
        ImGui::PopItemWidth();
        ImGui::EndDisabled();

        // Status Text
        ImGui::Text(rendering ? "Status: Rendering..." : "Status: Ready.");

        ImGui::End(); // End Controls window

        window.clear(sf::Color(20, 20, 20)); // Dark background for the window

        // --- Display the rendered image ---
        // sf::Vector2f and setPosition typically take floats.
        sf::RectangleShape image_background(sf::Vector2f(static_cast<float>(render_image_width), static_cast<float>(render_image_height)));
        image_background.setPosition(static_cast<float>(gui_panel_width), 0.0f); // Position next to GUI panel
        image_background.setFillColor(sf::Color(0, 0, 0)); // Black background for the render area
        window.draw(image_background);

        {
            std::lock_guard<std::mutex> lock(pixel_buffer_mutex); // Lock while accessing shared buffer/texture
            // Draw sprite only if pixel_buffer_rgba is not empty and texture is valid
            if (!pixel_buffer_rgba.empty() && rendered_texture.getSize().x > 0) {
                rendered_sprite.setPosition(static_cast<float>(gui_panel_width), 0.0f); // Position next to GUI panel
                window.draw(rendered_sprite);
            }
            if (rendering) {
                // Display "Rendering..." text overlay if in progress
                ImGui::SetNextWindowPos(ImVec2(static_cast<float>(gui_panel_width), 0.0f));
                ImGui::SetNextWindowSize(ImVec2(static_cast<float>(render_image_width), static_cast<float>(render_image_height)));
                ImGui::Begin("##RenderOverlay", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoInputs | ImGuiWindowFlags_NoBackground);
                // Center text manually
                ImGui::SetCursorPosX((static_cast<float>(render_image_width) - ImGui::CalcTextSize("Rendering...").x) / 2.0f);
                ImGui::SetCursorPosY((static_cast<float>(render_image_height) - ImGui::CalcTextSize("Rendering...").y) / 2.0f);
                ImGui::Text("Rendering...");
                ImGui::End();
            }
        }

        ImGui::SFML::Render(window); // Render ImGui elements
        window.display(); // Display everything
    }

    ImGui::SFML::Shutdown(); // Clean up ImGui-SFML
    return 0;
}-