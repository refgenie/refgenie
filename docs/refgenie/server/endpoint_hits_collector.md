# Endpoint Hits Collector in RefgenieServer

## Overview

The `endpoint_hits_collector` is a core component in RefgenieServer designed to track and process hits to specific API endpoints. Its primary use case is to collect download counts for assets, but it is built to be easily extensible for other analytics or logging purposes.

## How It Works

### Middleware Collection

- The `EndpointHitCounterMiddleware` is registered as a middleware in the FastAPI app.
- On each request, it inspects the matched route and path parameters.
- It records a hit by calling `record_hit` on the `endpoint_hits_collector` (attached to `app.state`).

### Background Processing

- The `endpoint_hits_collector` is periodically processed by a background job, managed by an APScheduler `BackgroundScheduler`.
- The scheduler is started in the FastAPI app's `lifespan` context and runs the `handle_endpoint_collector_hits` function at a configurable interval (default: 60 seconds)
- This function executes all registered handlers (see below), passing them the collected hits, and then resets the collector for the next interval.

### Current Usage: Download Counts

- The main use case in RefgenieServer is to collect and persist download counts for assets.
- Handlers in `refgenieserver.stats.handlers.HANDLERS_CATALOG` process the collected hits and update the relevant statistics (e.g., incrementing download counters in a database).

## Extending Functionality: Adding More Handlers

Developers can extend the analytics and processing capabilities by adding new handlers. Here’s how:

1. **Create a Handler Function**
   - A handler is a function that takes two arguments: the `EndpointCollector` instance and the `Refgenie` instance.
   - It should process the collected hits as needed (e.g., log to a file, update metrics, trigger alerts).

   ```python
   def my_custom_handler(endpoint_collector, refgenie):
       # Process collected hits
       hits = endpoint_collector.get_hits()
       # ... custom logic ...
   ```

2. **Register the Handler**
   - Add your handler to the `HANDLERS_CATALOG` dictionary in `refgenieserver/stats/handlers.py`:

   ```python
   HANDLERS_CATALOG = {
       "download_count": download_count_handler,
       "my_custom": my_custom_handler,
   }
   ```

3. **Automatic Execution**
   - All handlers in the catalog are automatically called by `handle_endpoint_collector_hits` at each interval.
   - No further changes are needed to the scheduler or main app.

## Configuration

- The interval for processing hits is controlled by the `DOWNLOAD_COUNT_DUMP_JOB_INTERVAL_SECONDS` environment variable (default: 60 seconds).
- This can be set in your environment or configuration files.

## Summary

- The `endpoint_hits_collector` provides a flexible, centralized way to collect and process endpoint usage data.
- It is managed automatically by middleware and a background scheduler.
- Developers can easily extend its functionality by adding new handlers to the handler catalog.
