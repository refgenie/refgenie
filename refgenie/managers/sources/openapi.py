from enum import Enum

from pydantic import BaseModel


class HttpMethod(Enum):
    GET = "GET"
    POST = "POST"
    PUT = "PUT"
    DELETE = "DELETE"
    HEAD = "HEAD"


class Tag(BaseModel):
    name: str
    description: str


class Endpoint(BaseModel):
    operationId: str


class EndpointsByMethod(BaseModel):
    get: Endpoint | None = None
    post: Endpoint | None = None
    put: Endpoint | None = None
    delete: Endpoint | None = None
    head: Endpoint | None = None


class Info(BaseModel):
    title: str
    version: str
    description: str


class RefgenieserverOpenApiSpec(BaseModel):
    openapi: str
    info: Info
    paths: dict[str, EndpointsByMethod]
    tags: list[Tag]
