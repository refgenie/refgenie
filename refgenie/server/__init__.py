"""
Production server module for refgenie.

This module provides the FastAPI application for running a production refgenie server
with analytics tracking, GA4GH DRS compliance, and data channel aggregation.

Usage:
    refgenie serve --port 8000
    # or
    uvicorn refgenie.server.main:app --host 0.0.0.0 --port 8000
"""
